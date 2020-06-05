      subroutine ini_st
!     --------------------
!     DESCRIPTION

!     Initializes density for model using an analytic 
!     density function, derived from an N^2 profile
!     Across-front use the TANH profile with tightness
!     as a measure of frontal spread
!     Larger the factor, tighter the front and larger b_xx.

!     ____________________
!     --------------------
!     SET VARIABLES

      USE header 
      implicit none 
      REAL(kind=rc_kind)  :: alpha, T_ocean, sbkgrnd,dtemp, &
           n2bkgrnd,dTdz,wiggles, amplitude,thy,grav,dz,dep_front
!      include 'header.f'
      parameter (alpha= 1.d-4, grav=10.)
      integer  i,j,k,n
!     yfront is the distance of the front from the coast
      yfront = 0.5*(yc(NJ+1) +yc(0))
!     ____________________
!     --------------------
!     DECLARE VARIABLES
!      T_ocean = 4.0 !T_ocean is the deep temperature    
      T_ocean = 6.0 !T_ocean is the surface temperature    
      sbkgrnd = 35.0 ! constant salinity
!      tightness = 0.03d0  
      tightness = 0.1d0  

      n2bkgrnd = 1.d-6    !N^2 value
      dTdz = n2bkgrnd/(alpha* grav)
      
      !set salnity constant,  front is only in temperature
      s = sbkgrnd
!-------------------------------------------------------------
      ! cross-front gradients
      drho = 0.02   !this is delta rho across the front
      dtemp=  drho/ ( alpha * 1000.)    ! this is delta temp , e.g 1 deg C
!     --------------------
      n=0
      !Set up constant stratification
!      T(:,:,0,n) = T_ocean
      T(:,:,NK+1,n) = T_ocean
      do j=0,NJ+1
         do i=0,NI+1
            do k=NK,0,-1
               dz = (zc(i,j,k+1) - zc(i,j,k))*DL
               T(i,j,k,n) = T(i,j,k+1,n) - dtdz* dz
            end do
         end do
      end do

!     --------------------
!     ADD FRONT
!     A frontal region is introduced to the channel.
 
      dep_front= 0.6*sum(D)/dble((NI+2)*(NJ+2))
      do j=0,NJ+1
         do i=0,NI+1
 
!     WIGGLE in FRONT
            wiggles=1.d0
            amplitude= 0.0d0
            yfront= 0.5d0*(yc(NJ+1) +yc(0))       & 
                 + amplitude*dsin(2.d0*PI*xc(i)/  & 
                 (0.5*(xc(NI+1)+xc(NI))/wiggles)  )
           thy = tanh(tightness*(yc(j)-yfront)*PI)
           do k=0,NK+1
!               z= DL*zc(i,j,k)
!                  if (z.ge.-z1) then 
!                     slfacnew = slfac
!                  else if (z.ge.-z2) then
!                     slfacnew = slfac*(z+z2)/(z2-z1)
!                  else 
!                     slfacnew = 0.d0
!                  end if
!                  s(i,j,k,0)= slfacnew*thy + s(i,NJ,k,0)
              if (zc(i,j,k).lt.dep_front) then
                 T(i,j,k,0) = ((-zc(i,j,k)+dep_front)/dep_front)*dtemp*thy + T(i,NJ,k,n)
              else
                 T(i,j,k,0) = T(i,NJ,k,n)
              end if
!               end if
            end do
         end do
      end do



      return
      end

