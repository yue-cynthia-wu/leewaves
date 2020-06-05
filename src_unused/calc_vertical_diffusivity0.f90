subroutine calc_vertical_diffusivity ! Calculates vertical diffusivity for tracers and momentum  KzTr and KzMom (global arrays)
  !     ---------------------------------------------                     
! Updates KzMom and KzTr (global variables)
#include "cppdefs.h"

  USE header
  !  Find kappa either from Gotm or using our own scheme.
  implicit none 
  !REAL(kind=rc_kind) :: dSin,posneg,facreverse
  integer i,j,k,m,step,nday
  integer iv_compute_kzl 
  REAL(kind=rc_kind) :: zdep,yw,yset,ymid,ustar,Ekdepth,ff,dissrate,bz,rdz,kapp,kappaz,edepinv,const
                                                                   
  REAL(kind=rc_kind), parameter :: Kzmax= 1.d-3, KzmaxTr=1.d-3 !  The viscosity is Kz*Kzmax 

  INTEGER OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS,CHUNK,NPROC

 !-------------------
 ! COMPUTATION OF KZ
 !-------------------
 ! Call windstress and heatflux before
 ! calc_vertical_diffusivity is called only once per timestep at the beginning of momentum

!#ifdef gotm_call
!    do j=1,NJ    
!    do i=1,NI       
!       call shearn2(i,j)    
!      call couple_gotm(i,j)       !  updates  KzMom and KzTr for the i,j column 
!    enddo
!    enddo
!#else
!  KzMom=1.d-5
!  KzTr = 1.d-5
!  return
! edepinv  = 0.04    !for run1 -  - decay scale of  25 m
  edepinv  = 0.08    !for run 2  - decay scale of 12.5 m 
  const= (-gpr*10.d0/R0)

  do j=1,NJ 
     do i=1,NI 

        ! Scheme 2
        !============
        do k=NK-1,1,-1
           zdep = zf(i,j,k)*DL
!run1         dissrate= (0.1d0 +stress_top(i,j)/0.1d0)* 10.d0**(0.04d0*zdep)*1d-7   ! m^2/s^3
!run2 make dropoff sharper
           dissrate= (0.1d0 +stress_top(i,j)/0.1d0)* 10.d0**(edepinv*zdep)*1d-7   ! m^2/s^3
           rdz =(rho(i,j,k+1) -rho(i,j,k))*wz(i,j,k)*DLinv
!          bz = dmax1(1d-10, (-gpr*10.d0/R0)*rdz )    ! run1 
           bz = dmax1(1d-6, const*rdz )     ! let the max bz be 1e-6
!          kapp= dmax1(0.2d0*dissrate/bz, 1.d-5)     ! m^2/s     ! run 1 
           kapp= dmax1(0.2d0*dissrate/bz, 5.d-7)     ! m^2/s      ! make min kappa=5.e-7 for BoB run2
           kappaz= dmin1(1.d-3,kapp)                 ! Kz Not to exceed 1e-3
!==     Tanh profile in Kz,  Use KzmaxTr=1.d-3                
!          thy = (1.d0 +tanh(((zdep-ztransit)/zextent)*PI))         &  
!     &             *0.5d0                                                           
!     thy varies from +1 to 0                                                 
!          KzMom(i,j,k)= kappaz/Kzmax   ! let's do away with Kzmax
           KzMom(i,j,k)= kappaz
           KzTr(i,j,k)= KzMom(i,j,k)
        end do
        KzMom(i,j,NK)= KzMom(i,j,NK-1)     ! value at k=NK not used, but filled for better plots
        KzTr(i,j,NK)= KzTr(i,j,NK-1)

!         ! Scheme 1
!         !===============
!          ustar = sqrt( stress_top(i,j) /R0 ) 
!          ff= ffc(i,j)*FPAR 
!          Ekdepth= 0.4d0*ustar/ff 
!          zextent= 0.5d0*Ekdepth !     Set zextent
!          ztransit = -Ekdepth 
!          do k=NK,0,-1 
!             Kz(i,j,k)= 1.d0 
!             zdep = zf(i,j,k)*DL 
!             thy = (1.d0 +tanh(((zdep-ztransit)/zextent)*PI))*0.5d0                                              
!             Kz(i,j,k)= max(0.01d0,thy)*KzmaxTr
!             KzMom(i,j,k)= Kz(i,j,k)
!             KzTr(i,j,k)= Kz(i,j,k)
!          end do
     enddo
  enddo
!#endif

     
return
end
