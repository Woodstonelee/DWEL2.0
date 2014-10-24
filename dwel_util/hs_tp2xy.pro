function hs_tp2xy,nxy,theta,phi,x,y
  compile_opt idl2
  
  ;this function returns the x,y coordinates
  ;from input theta,phi coordinates
  ;the inputs and outputs are assumed to be vectors
  
  status=0
  x=fltarr(nxy)
  y=fltarr(nxy)
  
  x=2.0*!dtor*theta/!pi
  
  y=x*cos(!dtor*phi)
  x=x*sin(!dtor*phi)
  
  return, status
  
end

