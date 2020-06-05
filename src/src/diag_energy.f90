subroutine diag_energy(step) 
!---------------------------------------------------                    
  USE header
                                                                       
  integer :: step
  REAL(kind=rc_kind) ::  const,enkin,umax,vmax,wmax,volmax
  logical isnan

  !sufficient to calculate every 10 time steps
  if (mod(step,10)==0) then
                                                                   
  const= EPS*EPS*delta*delta 

  enkin = SUM( Jac(1:NI,1:NJ,1:NK) * ( u(1:NI,1:NJ,1:NK,0)**2 + v(1:NI,1:NJ,1:NK,0)**2 + const*w(1:NI,1:NJ,1:NK,0)**2  ))
  volmax= maxval(Jac(1:NI,1:NJ,1:NK))
  umax= maxval(u(1:NI,1:NJ,1:NK,0))
  vmax= maxval(v(1:NI,1:NJ,1:NK,0))
  wmax= maxval(w(1:NI,1:NJ,1:NK,0))

  print*, "max u,v,w,Jac, #total kinetic energy = ", step, umax,vmax,wmax,volmax,enkin
  
  if (isnan(volmax)) stop
  if (isnan(umax)) stop
  if (isnan(vmax)) stop
  if (isnan(wmax)) stop
  if (isnan(enkin)) stop

  end if

  return 
END subroutine diag_energy                                       


logical function isnan(a) 
  USE header
  REAL(kind=rc_kind) ::a
  if (a.ne.a) then 
     isnan = .true. 
  else 
     isnan = .false. 
  end if
  return 
end 
