function dwel_eff_nu, wavelength, div, range, rtef
compile_opt idl2
;;
;; Calculate EVI telescope efficiency based on exponential model
;; Input parameters are:
;; div : EVI beam divergence (mrad)
;; range : an array of target ranges (metres) where efficiency is to be
;;           calculated
;; rtef is the range of the exponential
;; Returned value is an array with same dimensions as range.
;;
;; Set constants
;;
;mirror_diam = 100.0
;central_obs = 30.0
;focal_length = 100.0
;det_diam = 1.5
;
;;r_scale allows for the beam divergence in the range of the
;;telescope efficiency factor
;r_scale=1.0
;r_scale=(0.062338*(div^2.40472)+9.74752)/12.75
;
;;print,'r_scale=',r_scale
;r_max=max([mirror_diam,focal_length])/1000.0
;; Find valid ranges (> focal_length). Set result to 0 for invalid ranges.
;valid = where(range gt r_max, nvalid, compl=invalid, ncompl=n)
;eff_nu = replicate(1.0, n_elements(range))
;if (n gt 0) then eff_nu[invalid] = 0.0
;
;if (nvalid gt 0) then begin
;  ; Subset to valid ranges
;  rvalid = range[valid]
;
;  ; Check exponent
;  od = (rvalid/(rtef*r_scale))^2
;
;  ; Test relative size of image and od and apply efficiency relationship
;  ; as required.
;  w = where(od lt 30.0, nw)
;  if (nw gt 0) then eff_nu[valid[w]] = 1.0-exp(-od[w])
;  nan = where(finite(od) eq 0, nn)
;  if (nn gt 0) then eff_nu[valid[nan]] = 0.0
;  od=0b
;  rvalid = 0b
;endif
;eff_nu = ((eff_nu > 0.0) < 1.0)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;simplify DWEL telescope efficiency as usual one since we do not have 
;the efficience characteristic data now.
if (wavelength eq 1064) then begin
  c=5.43
endif
if (wavelength eq 1548) then begin
  c=5.87
endif

valid = where(range gt 0.01, nvalid, compl=invalid, ncompl=n)
eff_nu = replicate(1.0, n_elements(range))
if (n gt 0) then eff_nu[invalid] = 0.0
eff_nu[valid] = 1-exp(-(range[valid]/c)^2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

return, eff_nu
end
