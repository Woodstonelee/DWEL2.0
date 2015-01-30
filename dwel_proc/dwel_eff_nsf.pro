function dwel_eff_nsf, wavelength, range, par=par
  compile_opt idl2
  ;
  ;Calculate DWEL telescope efficiency based on growth model
  ;Input parameters are:
  ;wavelength : the DWEL band
  ;range : an array of target ranges (metres) where efficiency is to be
  ;           calculated
  ;keyword argument, par: calibration parameters
  ; Returned value is an array with same dimensions as range.
  ;
  ; k(r) model function in latex:
  ; P_r(r) = \frac{C_0 \cdot \rho}{r^b} \cdot \frac{1}{ \left(
  ; 1+C_1\cdot e^{-C_2\cdot(r+C_3)} \right)^ {C_4} }

  if n_elements(par) ne 0 or arg_present(par) then begin
    c1 = par[0]
    c2 = par[1]
    c3 = par[2]
    c4 = par[3]
  endif else begin
    if (wavelength eq 1064) then begin
      c1=3413.743d0
      c2=0.895d0
      c3=15.640d0
      c4=c1
    endif
    if (wavelength eq 1548) then begin
      c1=5.133d0
      c2=0.646d0
      c3=1.114d0
      c4=c1
    endif
  endelse 
  zero=1.0d0/(1.0d0+double(c1)*exp(-double(c2)*double(c3)))^c4
  num_range=n_elements(range)
  valid = where(range gt 0.05, nvalid, compl=invalid, ncompl=n)
  eff_nsf = cmreplicate(1.0,num_range)
  if (n gt 0) then eff_nsf[invalid] = float(zero)
  if (nvalid gt 0) then begin
    logK=dblarr(nvalid)
    logK=cmreplicate(double(c2),nvalid)*(double(range[valid])+cmreplicate(double(c3),nvalid))
    pos=where(logK gt 30.0d0,npos)
    if (npos gt 0) then logK[pos]=30.0d0
    logK=cmreplicate(double(c4),nvalid)*alog(cmreplicate(1.0d0,nvalid)+cmreplicate(double(c1),nvalid)*exp(-logK))
    eff_nsf[valid] = float(exp(-logK))
  endif
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  pos=where((eff_nsf lt 0.0) or (eff_nsf gt 1.0),npos)
  if (npos gt 0) then print,'out of range in dwel_eff!'
  return, eff_nsf

  ;; compile_opt idl2
  ;; ;
  ;; ;Calculate DWEL telescope efficiency based on EVI empirical model
  ;; ;Input parameters are:
  ;; ;wavelength : the DWEL band
  ;; ;range : an array of target ranges (metres) where efficiency is to be
  ;; ;           calculated
  ;; ; Returned value is an array with same dimensions as range.
  ;; ;
  ;; if (wavelength eq 1064) then begin
  ;;   ck = 19.6294
  ;; endif
  ;; if (wavelength eq 1548) then begin
  ;;   ck = 27.1794
  ;; endif

  ;; num_range=n_elements(range)
  ;; valid = where(range gt 0.05, nvalid, compl=invalid, ncompl=n)
  ;; eff_nsf = cmreplicate(1.0,num_range)

  ;; if (n gt 0) then eff_nsf[invalid] = 0.0

  ;; if (nvalid gt 0) then begin
  ;;   rvalid = range[valid]

  ;;   ;; check exponent
  ;;   od = rvalid^2/ck
  ;;   w = where(od lt 30.0, nw)
  ;;   if (nw gt 0) then eff_nsf[valid[w]] = 1.0 - exp(-od[w])
  ;;   nan = where(finite(od) eq 0, nn)
  ;;   if (nn gt 0) then eff_nsf[valid[nan]] = 0.0
  ;; endif
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; pos=where((eff_nsf lt 0.0) or (eff_nsf gt 1.0),npos)
  ;; if (npos gt 0) then print,'out of range in dwel_eff!'
  ;; eff_nsf = ((eff_nsf > 0.0) < 1.0)
  ;; return, eff_nsf
end
