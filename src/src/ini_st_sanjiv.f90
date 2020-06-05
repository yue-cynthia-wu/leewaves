      subroutine stprofile
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
!      implicit logical (a-z)
      implicit none 
!      include 'header.f'
      integer  i,j,k,trigger,m,TP
      double precision z, analyticeval, sbkgrnd, rhobkgrnd, tbkgrnd, T_ocean, z1, z2, slfac
      double precision thy, slfacnew, amplitude, widbar
      double precision nml, n21, n2dep, wiggles, mmm, ttt
      double precision  tt, qq, slope, yinter, alpha

      ! alpha is thermal expansion coeff to convert to temp   dtemp= drho /(rho0 *alpha)

!     yfront is the distance of the front from the coast
      yfront = 0.5*(yc(NJ+1) +yc(0))
!     ____________________
!     --------------------
!     DECLARE VARIABLES
      T_ocean = 4.0 !T_ocean is the ambient temperature    
      sbkgrnd = 35.0 ! constant salinity

      tightness = 0.03d0  

      mldepth = 105.d0     !Mixed layer depth
      nml = 1.0d-6         !N^2 value in mixed layer
      n21=  1.d-6   !n21 = 6.0d-5         !N^2 value at baroclinic depth
      n2dep = 1.d-6   != 1.5d-5       !N^2 value below depth of baroclinicity
      widbar = 10.d0       !Value used to determine width

!-----------when front is confined to the mixed layer----------
!       z1 = mldepth     ! Depth over which frontal strength is constant
!       z2 = mldepth +30 ! Max depth to which front extends
!------------when front extends below the mixed layer-----------

      z1 = mldepth       !Depth in which front is constant
      z2 = mldepth + 3.d0*mldepth  !Max depth to which front extends
!-------------------------------------------------------------
      slfac = 0.1   !this is delta rho
      trigger = 0  !0 means front trigged
!      wiggles = 10.d0
!      amplitude = 0.0d0
!     ____________________
!     --------------------
!     SET BACKGROUND DENSITY USING ANALYTIC PROFILE
      write(6,*) "tightness= ", tightness
      mmm = -zc(10,10,7)
      ttt = analyticeval(mldepth,mmm,nml,n21,n2dep,widbar)
      tt = mldepth + 13*widbar-1
!      qq = mldepth + 13*widbar


      do m=0, NK+1
      z=-zc(10,10,m)*1000
      TP = m-2
      if (z.le.tt) goto 100
      end do
100   continue

      write(6,*) "TP =  ", TP 

      !set salnity constant
      s = sbkgrnd

      do j=0,NJ+1
         do i=0,NI+1
            do k=TP,NK+1
                z= -zc(i,j,k)*1000
!                write(6,*) z
!                if (k.ge.TP) then
!                  write(6,*) z
!                 write(6,*) k
                  rhobkgrnd = analyticeval(mldepth,z,nml,n21,n2dep,widbar)
!                  write(6,*) sbkgrnd
                  tbkgrnd = T_ocean +  (rhobkgrnd -R0)/(alpha*R0)
                  T(i,j,k,0) =  tbkgrnd

!                else 

!                  slope = (s(1,1,TP-1,0)-s(1,1,TP-3,0))  & 
!                       / ((zc(1,1,TP-1)-zc(1,1,TP-3))*(-1000))


!                 sbkgrnd = slope * z - slope * zc(i,j,TP)*(-1000)  & 
!                + s(i,j,TP,0)
!                 s(i,j,k,0) = sbkgrnd
!               end if
            end do
         end do
      end do

      slope = (T(10,10,TP+1,0)-T(10,10,TP+2,0))  &  
             / ((zc(10,10,TP+1)-zc(10,10,TP+2))*(-1000))
      write(6,*) "slope = ", slope
      write(6,*) "zc1010TP= ", zc(10,10,TP)*1000
      write(6,*) "temp1010TP= ", T(10,10,TP,0)

      yinter =  slope * zc(10,10,TP+1)*(1000) + T(10,10,TP+1,0)      
      write(6,*) "yinter = ", yinter
      do j=0,NJ+1
         do i=0,NI+1
            do k=0,TP+1  
                z= -zc(i,j,k)*1000
                  tbkgrnd = slope * z + yinter

!- slope * zc(10,10,TP+1)*(1000) & 
!                 - s(10,10,TP+1,0)
                  T(i,j,k,0) = tbkgrnd
!                  write(6,*) sbkgrnd

            end do
         end do
      end do



!     ____________________
!     --------------------
!     ADD FRONT


!     A frontal region is introduced to the channel.
!     The front extends to a depth of z2. The front
!     is strongest to depth z1 and decreases in 
!     strength from z1 to z2

      if (trigger.eq.0) then
 
      do j=0,NJ+1
         do i=0,NI+1
 
!     WIGGLE in FRONT
            wiggles=1.d0
            amplitude= 0.1d0
            yfront= 0.5d0*(yc(NJ+1) +yc(0))       & 
                 + amplitude*dsin(2.d0*PI*xc(i)/  & 
                 (0.5*(xc(NI+1)+xc(NI))/wiggles)  )
           thy = tanh(tightness*(yc(j)-yfront)*PI)
           do k=1,NK
               z= DL*zc(i,j,k)
!              if (z.gt.-2*mldepth) then
                  if (z.ge.-z1) then 
                     slfacnew = slfac
                  else if (z.ge.-z2) then
                     slfacnew = slfac*(z+z2)/(z2-z1)
                  else 
                     slfacnew = 0.d0
                  end if
                  s(i,j,k,0)= slfacnew*(thy-25d-1) + s(i,NJ,k,0)
!               end if
            end do
         end do
      end do

      else 
        if (trigger.gt.0) goto 101

      end if

101   continue

!     ___________________
!     -------------------
!     SET BOUNDARY CONDITIONS

!     North-south boundaries set to solid walls.
!     East-west boundaries are periodic.

      do k=0,NK+1
         do i=1,NI
            s(i,0,k,0)= s(i,1,k,0)
            s(i,NJ+1,k,0)= s(i,NJ,k,0)
         end do
         do j=0,NJ+1
            s(0,j,k,0)= s(NI,j,k,0)
            s(NI+1,j,k,0)= s(1,j,k,0)
         end do
      end do
!     ____________________


      return
      end

