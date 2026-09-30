module interpolation
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |March 2025       |New module to manage Lagrangian tracer particles
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines which carry out interpolation
  !!


  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
#ifdef mp  
  use mpi_transfers
#endif  
  implicit none

  real(rkind) :: delta_fp
  integer(ikind) :: N_fp
  real(rkind),dimension(:),allocatable :: y_fp
  integer(ikind),dimension(:),allocatable :: proc_fp,nneighbour_fp

  real(rkind) :: delta_lp
  integer(ikind) :: N_lp
  real(rkind),dimension(:),allocatable :: y_lp,x_lp,u_lp,v_lp,cxx_lp,cxy_lp,cyy_lp
  integer(ikind),dimension(:),allocatable :: proc_lp,nneighbour_lp
  real(rkind),dimension(:,:),allocatable :: interp_w_ij


contains
!! ------------------------------------------------------------------------------------------------
  subroutine initialise_flux_points
     !! Generate a set of points on a vertical line at x=0, and find their nearest neighbours
     integer(ikind) :: i,j
     real(rkind),dimension(:),allocatable :: neighbour_dist_tmp
     real(rkind) :: dst_tmp
     
     !! How many points:     
     N_fp = floor((ymax-ymin)/smin_global)
     delta_fp = (ymax-ymin)/dble(N_fp)

     !! Calculate the positions of the flux-points
     allocate(y_fp(N_fp))
     do i=1,N_fp
        y_fp(i) = ymin + half*delta_fp + dble(i-1)*delta_fp
     end do
     
     !! Identify the processor and nearest neighbour node of each flux point
     allocate(neighbour_dist_tmp(N_fp));neighbour_dist_tmp = verylarge
     allocate(nneighbour_fp(N_fp));nneighbour_fp = -1
     do i=1,np
        if(abs(rp(i,1)).gt.s(i)) cycle !! Only loop over nodes near x=0
        
        !! Loop over all flux-points
        do j=1,N_fp
           !! Calculate i-j distance squared
           dst_tmp = rp(i,1)*rp(i,1) + (rp(i,2)-y_fp(j))**two

           !! If this is nearest SO FAR, and is within 2s, update list
           if(dst_tmp.lt.neighbour_dist_tmp(j).and.dst_tmp.le.s(i)*s(i)) then
              neighbour_dist_tmp(j) = dst_tmp
              nneighbour_fp(j) = i
           end if
 
        
        end do           
     
     end do

     !! Set the processor
     allocate(proc_fp(N_fp));proc_fp=-1
     do i=1,N_fp
        j=nneighbour_fp(i)
        if(j.ne.-1.and.j.le.npfb) then !! If it has a nearest neighbour, and that neighbour is <npfb, it belongs
           proc_fp(i) = iproc
        end if
!        if(proc_fp(i).ne.-1) then
!           write(6,*) iproc,i,y_fp(i),nneighbour_fp(i),neighbour_dist_tmp(i),proc_fp(i)
!        end if
     end do


     deallocate(neighbour_dist_tmp)

     return
  end subroutine initialise_flux_points
!! ------------------------------------------------------------------------------------------------  
  subroutine calculate_volumetric_flux(vol_flux,flux_length)
     !! This routine calculates the volumetric flux over the line x=0
     real(rkind),intent(out) :: vol_flux,flux_length
     integer(ikind) :: i,j,k,ii
     real(rkind) :: dudx,dudy,flux_this_segment
     
     !! Loop over all flux points
     vol_flux = zero;flux_length=zero
     do i=1,N_fp

        !! Check if this flux point belongs to this processor
        if(proc_fp(i).eq.iproc) then

           !! Add to the flux length
           flux_length = flux_length + delta_fp

           !! Identify the nearest particle
           ii = nneighbour_fp(i)
           
           !! Calculate the u-velocity gradients of particle ii
           dudx = zero;dudy = zero
           do k=1,ij_count(ii)
              j=ij_link(k,ii)     
              
              dudx = dudx + u(j)*ij_w_grad(1,k,ii)
              dudy = dudy + u(j)*ij_w_grad(2,k,ii)
           end do
           dudx = dudx - u(ii)*ij_w_grad_sum(1,ii)
           dudy = dudy - u(ii)*ij_w_grad_sum(2,ii)
     
           !! Calculate the contribution to the volumetric flux of this segment
           flux_this_segment = delta_fp*(u(ii) - rp(ii,1)*dudx + (y_fp(i)-rp(ii,2))*dudy)
           
           !! Add to the total flux
           vol_flux = vol_flux + flux_this_segment
     
        end if
     end do
     
     !! Reduce across processors if necessary
#ifdef mp  
     call global_reduce_sum(vol_flux)
     call global_reduce_sum(flux_length)
#endif

     !! Scale by the characteristic length-scale
     vol_flux = vol_flux*L_char
     flux_length = flux_length*L_char
     
     return
  end subroutine calculate_volumetric_flux
!! ------------------------------------------------------------------------------------------------        
  subroutine initialise_interp_line
     use rbfs
     use svdlib
     !! Read in a file containing a set of points for an interpolation line.
     integer(ikind) :: i,j,ii,k,i1
     real(rkind),dimension(:),allocatable :: neighbour_dist_tmp
     real(rkind) :: dst_tmp,ff1,rad,x,y
     real(rkind),dimension(10,10) :: amat
     real(rkind),dimension(10) :: wvec,xvec,psivec
     
     open(unit=314,file="interp_points.in",status="old")
     
     if(iproc.eq.0) open(unit=315,file="./data_out/lp_xy.out")
     
     !! How many points:     
     read(314,*) N_lp
     allocate(y_lp(N_lp),x_lp(N_lp),u_lp(N_lp),v_lp(N_lp),cxx_lp(N_lp),cxy_lp(N_lp),cyy_lp(N_lp))
     do i=1,N_lp
        read(314,*) x_lp(i),y_lp(i)
     end do
     u_lp=zero;v_lp=zero
     cxx_lp=zero;cxy_lp=zero;cyy_lp=zero

     !! Processor 0 open file for output and output positions...
     if(iproc.eq.0) then
        do i=1,N_lp
           write(315,*) x_lp(i),y_lp(i)        
        end do
        flush(315)
        close(315)
        open(unit=316,file="./data_out/lp_u.out") 
        open(unit=317,file="./data_out/lp_v.out")                 
        open(unit=318,file="./data_out/lp_cxx.out")
        open(unit=319,file="./data_out/lp_cxy.out")                         
        open(unit=320,file="./data_out/lp_cyy.out")                                                                  
     endif
     
     !! Identify the processor and nearest neighbour node of each flux point
     allocate(neighbour_dist_tmp(N_lp));neighbour_dist_tmp = verylarge
     allocate(nneighbour_lp(N_lp));nneighbour_lp = -1
     do i=1,np
        
        !! Loop over all flux-points
        do j=1,N_lp
           !! Calculate i-j distance squared
           dst_tmp = (rp(i,1)-x_lp(j))**two + (rp(i,2)-y_lp(j))**two

           !! If this is nearest SO FAR, and is within 2s, update list
           if(dst_tmp.lt.neighbour_dist_tmp(j).and.dst_tmp.le.s(i)*s(i)) then
              neighbour_dist_tmp(j) = dst_tmp
              nneighbour_lp(j) = i
           end if
         
        end do           
     
     end do

     !! Set the processor
     allocate(proc_lp(N_lp));proc_lp=-1
     do i=1,N_lp
        j=nneighbour_lp(i)
        if(j.ne.-1.and.j.le.npfb) then !! If it has a nearest neighbour, and that neighbour is <npfb, it belongs
           proc_lp(i) = iproc
        end if
     end do

    !! Build 4th order interpolants...
    allocate(interp_w_ij(nplink,N_lp));interp_w_ij = zero   
    do i=1,N_lp
    
       !! Check if this point belongs to this processor
       if(proc_lp(i).eq.iproc) then
              
       !! Loop over neighbouring "particles"
       amat = zero
       ii = nneighbour_lp(i)
       do k=1,ij_count(ii)
          j=ij_link(k,ii)
          
          !! Find relative position
          x = rp(j,1) - x_lp(i)
          y = rp(j,2) - y_lp(i)
                         
          !! Scale x,y for hermite
          rad = sqrt(x*x + y*y)/h(ii)
          ff1 = Wab(rad)
          x=x/h(ii);y=y/h(ii)
                    
          !! Build Vector of Taylor monomials
          xvec(1) = one
          xvec(2) = x
          xvec(3) = y
          xvec(4) = (1.0/2.0)*x*x
          xvec(5) = x*y
          xvec(6) = (1.0/2.0)*y*y
          xvec(7) = (1.0/6.0)*x*x*x
          xvec(8) = (1.0/2.0)*x*x*y
          xvec(9) = (1.0/2.0)*x*y*y
          xvec(10)= (1.0/6.0)*y*y*y
                    
          !! Build vector of Basis functions
          wvec(1) = ff1*one
          wvec(2) = ff1*Hermite1(x)
          wvec(3) = ff1*Hermite1(y)
          wvec(4) = ff1*Hermite2(x)
          wvec(5) = ff1*Hermite1(x)*Hermite1(y)
          wvec(6) = ff1*Hermite2(y)
          wvec(7) = ff1*Hermite3(x)
          wvec(8) = ff1*Hermite2(x)*Hermite1(y)
          wvec(9) = ff1*Hermite2(y)*Hermite1(x)
          wvec(10)= ff1*Hermite3(y)
          
          !! Build matrix
          do i1=1,10
             amat(i1,:) = amat(i1,:) + xvec(i1)*wvec(:)   !! Contribution to LHS for this interaction
          end do   
                             
       end do  ! End neighbour loop
    
       
       !! Solve linear system to get Psi
       psivec = zero;psivec(1) = one
       call svd_solve(amat,10,psivec)       
       
       !! Second loop over neighbours
       do k=1,ij_count(ii)
          j=ij_link(k,ii)      
       
          !! Find relative position
          x = rp(j,1) - x_lp(i)
          y = rp(j,2) - y_lp(i)
          
          !! Scale x,y for hermite
          rad = sqrt(x*x + y*y)/h(ii)
          ff1 = Wab(rad)
          x=x/h(ii);y=y/h(ii)           
                            
          !! Build vector of Basis functions
          wvec(1) = ff1*one
          wvec(2) = ff1*Hermite1(x)
          wvec(3) = ff1*Hermite1(y)
          wvec(4) = ff1*Hermite2(x)
          wvec(5) = ff1*Hermite1(x)*Hermite1(y)
          wvec(6) = ff1*Hermite2(y)
          wvec(7) = ff1*Hermite3(x)
          wvec(8) = ff1*Hermite2(x)*Hermite1(y)
          wvec(9) = ff1*Hermite2(y)*Hermite1(x)
          wvec(10)= ff1*Hermite3(y)         
         
          !! Store wij_interp
          interp_w_ij(k,i) = dot_product(psivec,wvec)

          !! Evaluate properties at grid point i
          !us(i) = us(i) + up(j,1)*wij          
          
       end do
    
     
       end if
    end do



     deallocate(neighbour_dist_tmp)
    
     return
  end subroutine initialise_interp_line
!! ------------------------------------------------------------------------------------------------        
  subroutine line_interpolate
     !! This routine calculates properties at the N_lp points.
     integer(ikind) :: i,j,k,ii
     
     !! Initialise interpolation point values to zero
     u_lp=-verylarge;v_lp=-verylarge
     cxx_lp=-verylarge;cxy_lp=-verylarge;cyy_lp=-verylarge
     
     !! Loop over all interpolation points    
     do i=1,N_lp

        !! Check if this point belongs to this processor
        if(proc_lp(i).eq.iproc) then

           !! Zero accumulators this processor only!
           u_lp(i) = zero;v_lp(i) = zero
           cxx_lp(i)=zero;cxy_lp(i)=zero;cyy_lp(i)=zero

           !! Identify the nearest particle
           ii = nneighbour_lp(i)
           
           !! Calculate the tr(c) gradients of particle ii
           do k=1,ij_count(ii)
              j=ij_link(k,ii)     
              
              u_lp(i) = u_lp(i) + (u(j))*interp_w_ij(k,i)
              v_lp(i) = v_lp(i) + (v(j))*interp_w_ij(k,i)              
              cxx_lp(i) = cxx_lp(i) + cxx(j)*interp_w_ij(k,i)
              cxy_lp(i) = cxy_lp(i) + cxy(j)*interp_w_ij(k,i)
              cyy_lp(i) = cyy_lp(i) + cyy(j)*interp_w_ij(k,i)                            
              
           end do             
        end if
     end do
     
     !! Reduce across processors if necessary
#ifdef mp  
     do i=1,N_lp
        call global_reduce_max(u_lp(i))
        call global_reduce_max(v_lp(i))
        call global_reduce_max(cxx_lp(i))
        call global_reduce_max(cxy_lp(i))                                                        
        call global_reduce_max(cyy_lp(i))   
     end do
#endif

     !! Write output
     if(iproc.eq.0) then
        write(316,*) u_lp(:)
        write(317,*) v_lp(:)
        flush(316)
        flush(317)
        write(318,*) cxx_lp(:)
        write(319,*) cxy_lp(:)
        write(320,*) cyy_lp(:)                
     endif
    
     return
  end subroutine line_interpolate
!! ------------------------------------------------------------------------------------------------        
!! ------------------------------------------------------------------------------------------------        
  function Hermite1(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: Hres
     Hres = z
  end function Hermite1
  function Hermite2(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: Hres
     Hres = z*z - 1.0d0
  end function Hermite2
  function Hermite3(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: Hres
     Hres = z*z*z - 3.0d0*z
  end function Hermite3
!! ------------------------------------------------------------------------------------------------                
end module interpolation
