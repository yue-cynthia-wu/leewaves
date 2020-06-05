subroutine tracersource(n,dtimel) 
  !     -----------------------                                           
  USE header

  !     tracer sources                                                    
  !                                                                       
  integer  i,j,k,n,m
  REAL(kind=rc_kind) :: dtimel,fac,tau,tauinv,vol,swrmax,swrmaxinv,swrd,rh,a1,b1,c1,d1,e1,f1
!  REAL(kind=rc_kind) :: trinit,facmortal,phyconvert 
  !                                                                       
  !     consump is the uptake (or if negative, then addition) of tracer   

  swrmax= 10.   ! used to non-dim the swr so the growth rate will be  (swrd/swrmax)* (1/tau)
  swrmaxinv= 1./swrmax
  !tau is the time scale for consump = 1day = 86400 seconds
  tauinv = 1.d0/(1.d0*86400.d0)    ! tauinv is in (per second)
  fac= dtimel/(1.d0*86400.d0*UL/LEN) ! dtime is non-dim time, TL=LEN/UL
  ! fac is non-dim
  ! J_A is 0.6. J_lambda1 is very small. J_lambda2 = 20 m - so only the second term matters. 
  !  swrd = swr(j)*( J_A*exp(zc(i,j,k)*DL/J_lambda1) + (1 - J_A)*exp(zc(i,j,k)*DL/J_lambda2)   )   
  ! At z=20,  swrd= 0.12*swr/swrmax (approx). At z=50 m, swrd=  0.032*swr/swrmax  If I choose swrmax = 10, 
  ! then for swr=100w/m2, swrd= 1.2 x (1 per day) at 20 m,  swrd= 0.32 x (1 per day) at 50m for 100 w/m2


  ! Tracer 1  is Nutrient without any consumption
  ! Tracer 2 is Nutrient with consumption
  ! Tracer 3 is Oxygen, reset to orig value at surface. 

  it=2
  do k=0,NK+1 
     do j=0,NJ+1 
        do i=0,NI+1 
           vol= Jac(i,j,k)*LEN*LEN*DL   !vol in m^3

           swrd = swr(j)*( J_A*exp(zc(i,j,k)*DL/J_lambda1) + (1. - J_A)*exp(zc(i,j,k)*DL/J_lambda2)   )*swrmaxinv
           if (Tr(it,i,j,k,n).gt.0.0) then
              ! consump in milimoles/ sec
              consump(i,j,k,1)= vol*tauinv*swrd*Tr(it,i,j,k,n)                           
!                (Tr(it,i,j,k,n)- trinit(j,k))  !Tr in milimoles per m^3
              Tr(it,i,j,k,n)= Tr(it,i,j,k,n) - fac*swrd*Tr(it,i,j,k,n)
           end if
        end do
     end do
  end do
  !
  ! prod is consump averaged in the x direction
  ! Multiply by 1.d-3 to get this in MOLES per SECOND
  do k=1,NK
     do j=1,NJ
        prod(j,k)= 0.d0
        do i=1,NI
           prod(j,k)= prod(j,k)+ consump(i,j,k,1)
        end do
        prod(j,k)= prod(j,k)*1.d-3/dble(NI)
     end do
  end do
  !                                                                       
! Tracer 3 Oxygen.  Reset to value at sig=21.5  at surface, 
  it= 3
  a1= -0.00404265970036948
  b1=  0.130618679680367
  c1= -0.914967639739012
  d1= 2.41698357406357
  rh= 0.
  k=NK
  Tr(it,:,:,NK,0)= a1*rh**3 + b1*rh**2 + c1*rh + d1




  do m=0,1 
     do k=0,NK+1 
        do j=0,NJ+1 
           do it=1,ntr 
              tr(it,0,j,k,m)= tr(it,NI,j,k,m) 
              tr(it,NI+1,j,k,m)= tr(it,1,j,k,m) 
           end do
        end do
     end do
  end do
  !                                                                       
  return 
end subroutine tracersource
