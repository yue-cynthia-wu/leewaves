subroutine tracerinit(stepl) 
  USE header  !only : Tr,rho
  INTEGER :: i,j,k,itr,stepl
  REAL(kind=rc_kind) :: a1,b1,c1,d1,rh

!     initializes tracer fields                                         
!     TRACER 1   NITRATE
!     =========                                                         
  a1= -0.08245
  b1= 5.072
  c1= -94.19
  d1= 499.972

!  a1= 0.318434421317965
!  b1= -3.38186354000513
!  c1= 10.8556165429202
!  d1= -3.36064810012217
!  e1= 2.32913778228182
!  rh= rho(i,j,k)-1021.5d0                    ! between rh = 0 and rh = 6

  do itr=1,2   ! Nitrate
  do k=1,NK
     do j=1,NJ
        do i=1,NI
           rh= rho(i,j,k)-1000.d0                    
           if (rh.gt.27.5) rh=27.5
           if (rh.le.21.5) rh= 21.5
!           if (rh.gt.6.) rh=6.
!           if (rh.le.0.) rh= 0.
!           Tr(itr,i,j,k,0)= a1*rh**4 + b1*rh**3 + c1*rh**2 + d1*rh + e1
           Tr(itr,i,j,k,0)=  a1*rh**3 + b1*rh**2 + c1*rh + d1

        end do
     end do
  end do
end do

  itr=3   ! Oxygen
!a =   -0.001307  (-0.001379, -0.001235)   From Mara  Feb 2016
!b =      0.1497  (0.1413, 0.1581)
!c =      -6.814  (-7.203, -6.424)
!d =         154  (145, 163)
!e =       -1730  (-1833, -1626)
!f =        7722  (7245, 8198)

  a1= -0.00404265970036948
  b1=  0.130618679680367
  c1= -0.914967639739012
  d1= 2.41698357406357


!  a1= -0.00130697707059244
!  b1= 0.149703248793906
!  c1= -6.81381630480589
!  d1= 154.032389314696
!  e1= -1729.58804695069
!  f1= 7721.53809567601

  do k=1,NK
     do j=1,NJ
        do i=1,NI
           rh= rho(i,j,k)-1021.5d0
           if (rh.gt.6.) rh=6.
           if (rh.le.0.) rh= 0.
!           if (rh.gt.27.5) rh=27.5
!           if (rh.le.21.5) rh= 21.5
           
!           Tr(itr,i,j,k,0)= a1*rh**5 + b1*rh**4 + c1*rh**3 + d1*rh**2 + e1*rh + f1
           Tr(itr,i,j,k,0)= a1*rh**3 + b1*rh**2 + c1*rh + d1

        end do
     end do
  end do

!      Tr = 0d0
!      Tr(itr,:,:,1,0) =1d0
end subroutine tracerinit
