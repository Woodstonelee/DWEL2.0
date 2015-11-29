function dwel_eff_oz, wavelength, range, par=par
  compile_opt idl2
  ;
  ;Calculate DWEL telescope efficiency based on growth model
  ;Input parameters are:
  ;wavelength : the DWEL band
  ;range : an array of target ranges (metres) where efficiency is to be
  ;           calculated
  ;keyword argument, par: parameters of telescope efficiency
  ; Returned value is an array with same dimensions as range.
  ;
  ; k(r) model function in latex:
  ; P_r(r) = \frac{C_0 \cdot \rho}{r^b} \cdot \frac{1}{ \left(
  ; 1+C_1\cdot e^{-C_2\cdot r} \right)^ {C_3} }

  if n_elements(par) ne 0 or arg_present(par) then begin
    c1 = par[0]
    c2 = par[1]
    c3 = par[2]
  endif else begin
    if (wavelength eq 1064) then begin
      c1=0.000071d0
      c2=0.383203d0
      c3=104574.246266d0
    endif
    if (wavelength eq 1548) then begin
      c1=0.000071d0
      c2=0.445864d0
      c3=104574.246266d0
    endif
  endelse
  zero=1.0d0/(1.0d0+double(c1)*exp(-double(c2)))^c3
  num_range=n_elements(range)
  valid = where(range gt 0.05, nvalid, compl=invalid, ncompl=n)
  eff_oz = cmreplicate(1.0,num_range)
  if (n gt 0) then eff_oz[invalid] = float(zero)
  if (nvalid gt 0) then begin
    logK=dblarr(nvalid)
    logK=cmreplicate(double(c2),nvalid) * (double(range[valid]))
    pos=where(logK gt 30.0d0,npos)
    if (npos gt 0) then logK[pos]=30.0d0
    logK=cmreplicate(double(c3),nvalid)*alog(cmreplicate(1.0d0,nvalid)+cmreplicate(double(c1),nvalid)*exp(-logK))
    eff_oz[valid] = float(exp(-logK))
  endif
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  pos=where((eff_oz lt 0.0) or (eff_oz gt 1.0),npos)
  if (npos gt 0) then print,'out of range in dwel_eff!'
  return, eff_oz
end
