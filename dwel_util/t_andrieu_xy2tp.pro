function t_Andrieu_xy2tp,nxy,x,y,theta,phi
compile_opt idl2

;this function returns the theta,phi coordinates
;from input x,y coordinates for the Andrieu projection
;the inputs and outputs are assumed to be vectors

common evi_constants,pi_val,rad,deg,eta

status=0
theta=fltarr(nxy)
phi=fltarr(nxy)

theta=x
phi=y

return, status

end

