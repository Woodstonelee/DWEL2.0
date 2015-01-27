function peak_int, x,y,x_int,y_int,offset
  compile_opt idl2
  
  ;this routine is rather specific to interpolating the peak to get range
  ;Basically, it assumes x and y have three elements with x as range and
  ;y as intensity of some kind. It assumes y[1]>= the other y values.
  ;In addition we assume x[0]<x[1]<x[2] and x_int is in the x range
  ;The peak needs the coefficient of x^2 at y[1] to be negative and
  ;we believe the interpolated value should be larger than the current maximum

  x = reform(x)
  y = reform(y)
  
  x_int=x[1]
  y_int=y[1]
  offset=0s
  
  ;solve for the quadratic pn
  d1=(y[2]-y[1])/(x[2]-x[1])
  d0=(y[1]-y[0])/(x[1]-x[0])
  c=(d1-d0)/(x[2]-x[0])
  
  if (c gt -1.0e-7) then return,1b
  
  b=d0-c*(x[1]+x[0])
  x_int=-b/(2.0*c)
  
  if (x_int lt x[0]-1.0e-7) then x_int=x[0]
  if (x_int gt x[2]+1.0e-7) then x_int=x[2]
  
  value=min(abs(x-x_int),offset)
  offset=offset-1s
  
  a=y[1]-b*x[1]-c*x[1]^2
  y_int=a+b*x_int+c*x_int^2
  
  return,0b
end

