subroutine mixing_vertical(var,vardif,m,step,iv_compute_kzl) 
  !     ---------------------------------------------                     
#include "cppdefs.h"

  USE header
  !     Wind Stress can be specified here                                 
  !     Kz*du/dz = -Tx/rho,   Kz*dv/dz = -Ty/rho                          
  !     use level m                                                       
  !     computes d2s/dz2 at the cell centers.                             
  !     fricu and fricv are u,v diffusion terms (in header)               
  !     that need to be saved (if m.eq.0) for n2budget)                   

  implicit none 
  !REAL(kind=rc_kind) :: dSin,posneg,facreverse
  integer i,j,k,m,step,nday
  integer iv_compute_kzl 
  REAL(kind=rc_kind) :: zdep,yw,yset,ymid,ustar,Ekdepth,ff                  
                                                                   
  REAL(kind=rc_kind), dimension(0:NI+1,0:NJ+1, 0:NK+1) :: var 
  REAL(kind=rc_kind) :: vardif(NI,NJ,NK)

  REAL(kind=rc_kind) :: dvardzfc(0:NK)
  REAL(kind=rc_kind) :: Kdudzb,Kdvdzb,Kdudzt,Kdvdzt,wzkth,fact,Cp     
  REAL(kind=rc_kind) :: fac,facb,dfactor,facy,rhoinv,diff
  REAL(kind=rc_kind) :: wgt,day,ztransit,zextent,thy                                         
  REAL(kind=rc_kind), parameter :: Kzmax= 1.d-3, KzmaxTr=1.d-3 !  The viscosity is Kz*Kzmax 

  INTEGER OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS,CHUNK,NPROC

  facb = RR*DL      ! Linear Drag                                                       
  ! facb = RR*DL*UL ! Quadratic drag                                                    

  fac= 1.d0/(UL*DL*delta) 
  fact = DL/UL 

 !-------------------
 ! COMPUTATION OF KZ
 !-------------------

#ifdef gotm_call
   do j=1,NJ    
      do i=1,NI       
         if (ivb==1) then
            call shearn2(i,j)    
         endif
         call couple_gotm(i,j)    

         if (selvar==1 .OR. selvar==2) then
            Kz(i,j,1:NK) = KzTr(i,j,1:NK)  
         else if (selvar==3 .OR. selvar==4 .OR. selvar==5) then
            Kz(i,j,1:NK) = KzMom(i,j,1:NK)  
         end if
       enddo
   enddo
#else
   if(iv_compute_kzl==1) then
      do j=1,NJ 
      do i=1,NI 
         ustar = sqrt( stress_top(i,j) /R0 ) 
         ff= ffc(i,j)*FPAR 
         Ekdepth= 0.4d0*ustar/ff 
         zextent= 0.5d0*Ekdepth !     Set zextent
         ztransit = -Ekdepth 
         do k=NK,0,-1 
            Kz(i,j,k)= 1.d0 
            zdep = zf(i,j,k)*DL 
            thy = (1.d0 +tanh(((zdep-ztransit)/zextent)*PI))*0.5d0                                              
            Kz(i,j,k)= max(0.01d0,thy)*KzmaxTr
         end do
      enddo
      enddo     
   else 
      if (selvar==1 .OR. selvar==2) then
      Kz(1:NI,1:NJ,1:NK) = KzTr
   else if (selvar==3 .OR. selvar==4 .OR. selvar==5) then
      Kz(1:NI,1:NJ,1:NK) = KzMom
   endif
   endif
#endif

!$OMP PARALLEL DO PRIVATE(j,i,k,dvardzfc,dfactor) SHARED(vardif) COLLAPSE(2)

  do j=1,NJ 
    do i=1,NI 

      !-------------------
      ! ADDITIONAL COMPUTATIONS
      !-------------------
                                                                      
      !     Quadratic drag                                                    
      !=            Kdudzb= facb*u(i,j,1,m)*abs(u(i,j,1,m))                  
      !=            Kdvdzb= facb*v(i,j,1,m)*abs(v(i,j,1,m))                  
      !     Linear Drag                                                       
      !Kdudzb= facb*u(i,j,1,m) 
      !Kdvdzb= facb*v(i,j,1,m) 
      !                                                              
      !rhoinv = 1.d0/rho(i,j,NK) 
      !Kdudzt= stressx(j)*rhoinv*fact 
      !Kdvdzt= stressy*rhoinv*fact 


      !-------------------
      ! COMPUTATION OF THE VARIATIONS OF THE VARIABLE ALONG Z
      !-------------------

      do k=1,NK-1 
        dvardzfc(k)= wz(i,j,k)*(var(i,j,k+1)-var(i,j,k)) 
      end do 

      dvardzfc(0)= 0.d0 
      dvardzfc(NK)= 0.0 


      !-------------------
      ! COMPUTATION OF THE FLUX DIVERGENCE
      !-------------------

      k=1 
!     ---
      dfactor=  fac*Jac(i,j,k)*wz(i,j,k) 
#ifdef implicit
      mat_B(i,j,k,selvar) = (-1.d0)*dfactor*(Kz(i,j,k)*wz(i,j,k) ) 
      mat_C(i,j,k,selvar) = dfactor*Kz(i,j,k)*wz(i,j,k)
#else
      vardif(i,j,k)= dfactor*(dvardzfc(k)*Kz(i,j,k)- dvardzfc(k-1)*Kz(i,j,k-1) )     
#endif

      do k=2,NK-1 
!     -----------
        dfactor=  fac*Jac(i,j,k)*wz(i,j,k) 

#ifdef implicit
	mat_A(i,j,k,selvar) = dfactor*kz(i,j,k-1)*wz(i,j,k-1)
	mat_B(i,j,k,selvar) = (-1.d0)*dfactor*( Kz(i,j,k)*wz(i,j,k) + Kz(i,j,k-1)*wz(i,j,k-1) )
	mat_C(i,j,k,selvar) = dfactor*Kz(i,j,k)*wz(i,j,k)
#else
        vardif(i,j,k)= dfactor*(dvardzfc(k)*Kz(i,j,k)- dvardzfc(k-1)*Kz(i,j,k-1) )    
#endif

      end do 
                                                                       
      k=NK 
!     ----                                                        
      dfactor=  fac*Jac(i,j,k)*wz(i,j,k) 

#ifdef implicit
      mat_A(i,j,k,selvar) = dfactor*Kz(i,j,k-1)*wz(i,j,k-1)
      mat_B(i,j,k,selvar) = (-1.d0)*dfactor*Kz(i,j,k-1)*wz(i,j,k-1)    
#else
      vardif(i,j,k)= dfactor*(dvardzfc(k)*Kz(i,j,k)- dvardzfc(k-1)*Kz(i,j,k-1) )       
#endif

    end do ! i
  end do ! j
!$OMP END PARALLEL DO                                                                   

return 
END                                           
                                                                        