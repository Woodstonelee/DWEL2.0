function peak_int_arr, x,y,x_int,y_int,offset
;+
;this is the matrix version of peak_int function.
;this routine is rather specific to interpolating the peak to get range
;Basically, it assumes x and y have three elements with x as range and
;y as intensity of some kind. It assumes y[1]>= the other y values.
;In addition we assume x[0]<x[1]<x[2] and x_int is in the x range
;The peak needs the coefficient of x^2 at y[1] to be negative and
;we believe the interpolated value should be larger than the current maximum
;;
;INPUTS:
;; x = array[3, n], three columns are x coordinates of three points for
;; interpolation. n rows are number of interpolations
;;
;; y = array[3, n], three columns are y coordinates of three points for
;; interpolation. n rows are number of interpolations
;;
;OUTPUTS:
;; x_int = vector, length of n, n interpolated x coordinates.
;;
;; y_int = vector, length of n, n interpolated y coordinates.
;;
;; offset = vector, length of n, the shift of x coordinate between interpolated
;; position and the closest bin. 
;;
;KEYWORDS:
;;
;RETURN:
;;
;MODIFICATION HISTORIES:
;-
;

  compile_opt idl2
  
  x = float(x)
  y = float(y)

  x_int=x[1, *]
  y_int=y[1, *]

  n = n_elements(x)
  if n ne n_elements(y) then begin
    print, 'peak_int_arr, x and y do not have the same number of interpolations!'
    return, 1b
  endif 

  offset=intarr(n)
  
  ;solve for the quadratic pn
  d1=(y[2, *]-y[1, *])/(x[2, *]-x[1, *])
  d0=(y[1, *]-y[0, *])/(x[1, *]-x[0, *])
  c=(d1-d0)/(x[2, *]-x[0, *])

  b = fltarr(n)
  a =  fltarr(n)

  badint_ind = where(c gt -1.0e-7, count, complement=goodint_ind, ncomplement=goodcount)
  if goodcount eq 0 then begin
    print, 'peak_int_arr, all are bad interpolations'
    return, 1b
  endif 
  
  b[goodint_ind]=d0[goodint_ind]-c[goodint_ind]*(x[1, goodint_ind]+x[0, goodint_ind])
  x_int[goodint_ind]=-b[goodint_ind]/(2.0*c[goodint_ind])
  
  tmpind = where(x_int lt x[0, *]-1.0e-7, tmpcount)
  if tmpcount gt 0 then begin
    x_int[tmpind] = x[0, tmpind]
  endif 

  tmpind = where(x_int gt x[2, *]+1.0e-7, tmpcount)
  if tmpcount gt 0 then begin
    x_int[tmpind] = x[2, tmpind]
  endif 
  
  value=min(abs(x[*, goodint_ind]-transpose(cmreplicate(x_int[goodint_ind], 3))),tmpminind,dimension=1)
  offset = array_indices(x[*, goodint_ind], tmpminind)
  offset = offset[0, *]
  offset=offset-1s
  
  a[goodint_ind]=y[1, goodint_ind]-b[goodint_ind]*x[1, goodint_ind]-c[goodint_ind]*x[1, goodint_ind]^2
  y_int[goodint_ind]=a[goodint_ind]+b[goodint_ind]*x_int[goodint_ind]+c[goodint_ind]*x_int[goodint_ind]^2

  if count gt 0 then begin
    x_int[badint_ind] = !values.f_nan
    y_int[badint_ind] = !values.f_nan
  endif 
  
  return,0b
end

