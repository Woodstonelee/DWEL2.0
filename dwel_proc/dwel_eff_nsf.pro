function dwel_eff_nsf, wavelength, range
compile_opt idl2
;
;Calculate DWEL telescope efficiency based on EVI empirical model
;Input parameters are:
;wavelength : the DWEL band
;range : an array of target ranges (metres) where efficiency is to be
;           calculated
; Returned value is an array with same dimensions as range.
;
if (wavelength eq 1064) then begin
  ck = 19.6294
endif
if (wavelength eq 1548) then begin
  ck = 27.1794
endif

num_range=n_elements(range)
valid = where(range gt 0.05, nvalid, compl=invalid, ncompl=n)
eff_nsf = cmreplicate(1.0,num_range)

if (n gt 0) then eff_nsf[invalid] = 0.0

if (nvalid gt 0) then begin
  rvalid = range[valid]
  
  ;; check exponent
  od = rvalid^2/ck
  w = where(od lt 30.0, nw)
  if (nw gt 0) then eff_nsf[valid[w]] = 1.0 - exp(-od[w])
  nan = where(finite(od) eq 0, nn)
  if (nn gt 0) then eff_nsf[valid[nan]] = 0.0
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pos=where((eff_nsf lt 0.0) or (eff_nsf gt 1.0),npos)
if (npos gt 0) then print,'out of range in dwel_eff!'
eff_nsf = ((eff_nsf > 0.0) < 1.0)
return, eff_nsf
end
