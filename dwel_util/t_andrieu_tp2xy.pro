function t_Andrieu_tp2xy,nxy,theta,phi,x,y
  compile_opt idl2
  
  ;this function returns the x,y coordinates for the Andrieu projection
  ;from input theta,phi coordinates. the inputs and outputs are assumed
  ;to be vectors
  
  common evi_constants,pi_val,rad,deg,eta
  
  status=0
  x=fltarr(nxy)
  y=fltarr(nxy)
  
  x = theta
  y = phi
  
  return, status
  
end
