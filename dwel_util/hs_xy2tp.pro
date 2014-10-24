function xy2tp,nxy,x,y,theta,phi
  compile_opt idl2
  
  ;this function returns the theta,phi coordinates
  ;from input x,y coordinates
  ;the inputs and outputs are assumed to be vectors
  
  status=0
  theta=fltarr(nxy)
  phi=fltarr(nxy)
  
  theta=sqrt(x^2+y^2)
  phi=atan(x,y)
  
  pos=where(phi lt 0.0,count)
  if (count gt 0) then begin
    phi=phi+2.0*!pi
  endif
  
  theta=theta*!pi/2.0
  theta=theta*!radeg
  phi=phi*!radeg
  
  return, status
  
end

