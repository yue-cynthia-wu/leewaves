   subroutine viscosity
!  ---------------------------------------------                     
   USE header
!  Compute a Richardson number based vertical viscosity              
!  coefficient Kz using the final velocities u,v and density         
!  with the n=0 subscript                                            
!  Kz is situated at k cell faces                                    
!  The algorithm is from Rev Geophys., p 373, 1994. Large et al.                                                                 
!  Assumes that rho is evaluated just prior to this call             

   integer i,j,k 
   integer, parameter :: (n1=3) 
   REAL(kind=rc_kind), parameter :: (RiCr= 0.7d0) 
   REAL(kind=rc_kind) :: grho,Ri
   REAL(kind=rc_kind), dimension(0:NI+1,0:NJ+1,0:NK+1) :: dudz,dvdz,drdz,vshear,bvfreq

   grho= 9.81/R0      ! g/rho_0
   DLinv = 1.0/DL 
	  
!  Set kz(k=0)= 1. It is required only for the s,T equations,
!  since for the momentum equations,  K*du/dz= ru.                         
!  Kz is not needed at k=0 and NK since stress bcs are used.         
   Kz(:,:, 0)= 0.d0  ! size of Kz is (NI,NJ,0:NK)
   Kz(:,:,NK)= 0.d0 
   
   ! ===========  Initialize  =========== !
!    dudz(:,:,:)  =0.d0
!    dvdz(:,:,:)  =0.d0
!    drdz(:,:,:)  =0.d0
!    vshear(:,:,:)=0.d0
!    bvfreq(:,:,:)=0.d0
   
   do k=1,NK-1 
      ! borrowed from "diag_vort.f90"
      dudz(:,:,k)= 0.5*( (  u(:,:,k+1,0)-  u(:,:,k-1,0))*wz(:,:,k) )
      dvdz(:,:,k)= 0.5*( (  v(:,:,k+1,0)-  v(:,:,k-1,0))*wz(:,:,k) )
      drdz(:,:,k)= 0.5*( (rho(:,:,k+1)  -rho(:,:,k-1)  )*wz(:,:,k) )
   end do
	 
   vshear(:,:,:)= ( dudz(:,:,:)*dudz(:,:,:)+dvdz(:,:,:)*dvdz(:,:,:) )*UL*UL/(DL*DL)		  
   bvfreq(:,:,:)= -grho*drdz(:,:,:)*DLinv ! N^2=-g/rho_0*(drho/dz) in s^-2 re-dim by DLinv 

   do i=1,NI
      do j=1,NJ
		  do k=1,NK-1
            if (vshear(i,j,k).eq.0.d0) then ! if (vshear.le.1.d-12) then                                    
               Ri= 100.d0 
            else 
               Ri= bvfreq(i,j,k)/vshear(i,j,k)
            end if 

            if (Ri.le.0.d0) then    ! unstable density profile => Ri < 0
               Kz(i,j,k)= 1.d0
            else if (Ri.lt.RiCr) then 
               Kz(i,j,k)= (1.d0 - (Ri*Ri)/(RiCr*RiCr))**n1 
            else 
               Kz(i,j,k)= 0.d0 
            endif 

         end do                                            
      end do 
   end do
   
!     The value of Kz to be used is  Kz(k) *Kzmax                       

   return 
   END                                           
