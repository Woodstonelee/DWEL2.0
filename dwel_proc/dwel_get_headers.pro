function DWEL_get_headers, in_fid, headers
compile_opt idl2
; routine to read all user-defined header info from EVI files.
; calling routine must have defined the headers structure with the
; tag 'f_base' containing the base name of the EVI file


; Initialise the headers_present counter
headers=create_struct('headers_present',0s,headers)

if (in_fid lt 0L) then return,0b

; Create lookup table of possible header values and a 'short name' used
; as a tag for keeping track of whether the values are present.
; These header info names should not be defined in the structure before
; calling this routine.
; The convention with the old version of this routine was to use the
; 'short name' with '_present' appended to log the presence of the given
; header record. This facility has been retained and and entry is required
; in this lookup table to register the short name for any new header records.
; This code should be identical to that in put_headers.pro
; The code below (between the ;*** marks) should be identical to that in
; put_headers.pro and is the only place where editing is required to add
; new header records.
; 
; Added three forgotten (so unwritten) header records 6/1/2011
; 
;***
header_info = { $
DWEL_scan_info:'run', $
DWEL_diagnostic_info:'diag', $
DWEL_base_fix_info:'base', $
DWEL_sat_info:'sat', $
DWEL_pfilter_info:'pfilter', $
DWEL_apprefl_info:'apprefl', $
DWEL_projection_info:'proj', $
DWEL_andrieu_zenith:'andzen',$
DWEL_andrieu_azimuth:'andaz',$
DWEL_cylindrical_projection_info:'cylproj', $
DWEL_cylindrical_rvals:'cylrvals', $
DWEL_cylindrical_azumith:'cylaz', $
DWEL_gap_dist_info:'gap', $
DWEL_zen_azm_average_info:'web', $
DWEL_trunk_info:'trunk',$
DWEL_pointcloud_info:'ptcld', $
DWEL_filtered_fix_info:'filtfix',$
dwel_adaptation:'dwel_import' $
}
;*** End of repeated code.

nnames = n_tags(header_info)
header_names = tag_names(header_info)

;get the DWEL run header information
for i=0, nnames-1 do begin
  info = envi_get_header_value(in_fid, header_names[i],undefined=undef)
  if (undef) then begin
      info='NO DWEL '+header_info.(i)+' information in header of '$
            +strtrim(headers.f_base,2)
      headers=create_struct(header_info.(i)+'_present',0b, headers)
  endif else begin
      headers=create_struct(header_info.(i)+'_present',1b, headers)
      num=DWEL_header_parse(info)
      headers.headers_present=headers.headers_present+1s
  endelse
  headers=create_struct(header_names[i],info,headers)
  info=0b
endfor

header_info=0b

return,1b

end
