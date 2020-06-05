    subroutine write_cdf_2D_y(jslice,counter_2d,n)
  !     ------------------------------------------------------------      
      USE header,ONLY : NI,NJ,NK,ntr,nconsume,s,T,rho,u,v,w,vor,xc,zc,DL, &
           time_seconds,time_days,dirout,rc_kind
  !     reads in a netcdf file and interpolates the data from the sigma gr
  !     onto a z level                                                    

#include "netcdf.inc" 
      integer :: itr,jslice,counter_2d,i,j,k,n
      REAL(kind=4) ::  xslice(NI),zslice(NI,NK),                &
     &     Tslice(NI,NK),rslice(NI,NK),uslice(NI,NK),   &
     &     vslice(NI,NK),wslice(NI,NK),vorslice(NI,NK)      
      REAL(kind=rc_kind) ::  rcode
                                                                        
      character(LEN=150) outname 

      integer start(4),count(4),counttim,dims(4),dims4(4),dimstim,       &
           start2d(2),count4(4),start4(4) 
      integer :: iddatfile,idvx,idvz,idvtim,idvu,idvv,idvw,idvvor,      &
           idvbz,idvby,idvrho,idvt,idvs 
      j=jslice                                    
      do k=1,NK 
         do i=1,NI 
            zslice(i,k)= zc(i,j,k)*DL 
            xslice(i)= xc(i) 
            !for tracer
!            do itr=1,ntr 
!               trslice(itr,i,k)= Tr(itr,i,j,k,n) 
!            end do 
!             Tslice(i,k)= T(i,j,k,n) 
            rslice(i,k)= rho(i,j,k) 
            uslice(i,k)= u(i,j,k,n) 
            vslice(i,k)= v(i,j,k,n) 
            wslice(i,k)= w(i,j,k,n) 
            vorslice(i,k)= vor(i,j,k) 
         end do 
      end do 

!     write yslice to a netcdf file      

      WRITE(outname,'("yslice_",I3.3,".cdf")') jslice                                      
      if (counter_2d.eq.1) then 
         idDatFile =  nccre(TRIM(dirout)//outname,NCCLOB,rcode) 

         dims(1) = ncddef(idDatFile,'x',NI,rcode) 
         dims(2) = ncddef(idDatFile,'y',1,rcode) 
         dims(3) = ncddef(idDatFile,'z',NK,rcode) 
         dims(4) = ncddef(idDatFile,'time',NCUNLIM,rcode) 

         dimstim= dims(4)

         !this is for tracer
         dims4(1) = dims(1) 
         dims4(2) = ncddef(idDatFile,'ntr',ntr,rcode) 
         dims4(3) = dims(3) 
         dims4(4) = dims(4) 

         idvx = ncvdef(idDatFile,'xc',NCFLOAT,1,dims(1),rcode) 
         idvz = ncvdef(idDatFile,'zc',NCFLOAT,3,dims,rcode) 
         idvtim = ncvdef(idDatFile,'day',NCFLOAT,1,dimstim,rcode) 
         idvt = ncvdef(idDatFile,'temp',NCFLOAT,4,dims,rcode)
!         idvt = ncvdef(idDatFile,'tr',NCFLOAT,4,dims4,rcode) 
         idvrho = ncvdef(idDatFile,'rho',NCFLOAT,4,dims,rcode) 
         idvu = ncvdef(idDatFile,'u',NCFLOAT,4,dims,rcode) 
         idvv = ncvdef(idDatFile,'v',NCFLOAT,4,dims,rcode) 
         idvw = ncvdef(idDatFile,'w',NCFLOAT,4,dims,rcode) 
         idvvor = ncvdef(idDatFile,'vor',NCFLOAT,4,dims,rcode) 

         CALL ncendf(idDatFile,rcode) 

      else 
         idDatFile = ncopn(TRIM(dirout)//outname, NCWRITE,rcode) 
         idvtim  = NCVID(idDatFile, 'day',RCODE) 
         idvt = NCVID(idDatFile, 'temp', RCODE) 
         idvrho = NCVID(idDatFile, 'rho', RCODE) 
         idvu = NCVID(idDatFile, 'u', RCODE) 
         idvv = NCVID(idDatFile, 'v', RCODE) 
         idvw = NCVID(idDatFile, 'w', RCODE) 
         idvvor = NCVID(idDatFile, 'vor', RCODE) 
      endif 

      write(6,*) 'in writecdf_yslice', 'counter,days', counter_2d,time_days
      count(1)= NI 
      count(2)= 1
      count(3)= NK
      count(4)= 1

      !this is for tracers
      count4(1)= NI 
      count4(2)= ntr 
      count4(3)= NK
      count4(4)= 1 

      counttim= 1

      start2d(1)= 1 
      start2d(2)= 1 

      start(1)= 1 
      start(2)= 1 
      start(3)= 1
      start(4)= counter_2d

      start4(1)= 1 
      start4(2)= 1 
      start4(3)= 1 
      start4(4)= counter_2d 

      if (counter_2d.eq.1) then 
         CALL ncvpt(idDatFile,idvx,start(1),count(1),xslice, rcode) 
         CALL ncvpt(idDatFile,idvz, start, count, zslice, rcode) 
      endif 
      CALL ncvpt(idDatFile,idvtim,start(4),counttim,time_days,rcode) 
!       CALL ncvpt(idDatFile,idvt, start, count, Tslice, rcode) 
      CALL ncvpt(idDatFile,idvrho, start, count, rslice, rcode) 
      CALL ncvpt(idDatFile,idvu, start, count, uslice, rcode) 
      CALL ncvpt(idDatFile,idvv, start, count, vslice, rcode) 
      CALL ncvpt(idDatFile,idvw, start, count, wslice, rcode) 
      CALL ncvpt(idDatFile,idvvor, start, count, vorslice, rcode) 

      CALL ncclos(idDatFile,rcode) 

      END                                           
