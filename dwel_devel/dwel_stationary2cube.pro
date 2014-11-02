;+
; NAME:
;DWEL_STATIONARY2CUBE
;
;
; PURPOSE:
;Fake a data cube from a stationary scan for DWEL program to process easily
;later. 
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;
;
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;-
pro dwel_stationary2cube, instathdf5file, outcubefile, nsample=nsample
  compile_opt idl2

  ;; check the keyword argument
  if n_elements(nsample) ne 0 or arg_present(nsample) then begin
    ;; user gives a nsample
  endif else begin
    ;; no nsample given, use default values
    nsample = 3142
  endelse 

  
end
