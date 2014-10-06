function DWEL_put_headers, in_fid, headers
compile_opt idl2

if (in_fid lt 0L) then return,0b

If (headers.headers_present le 0s) then return,1b

; Create lookup table of possible header values and a 'short name' used
; as a tag for keeping track of whether the values are present.
; These header info names should not be defined in the structure before
; calling this routine.
; The convention with the old version of this routine was to use the
; 'short name' with '_present' appended to log the presence of the given
; header record. This facility has been retained and and entry is required
; in this lookup table to register the short name for any new header records.
; The code below (between the ;*** marks) should be identical to that in
; get_headers.pro and is the only place where editing is required to add
; new header records.
; 
; Added two forgotten (so unwritten) header records 6/1/2011
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
DWEL_filtered_fix_info:'filtfix',$
DWEL_pointcloud_info:'ptcld', $
DWEL_adaptation:'dwel_import' $
}
;*** End of repeated code.

nnames = n_tags(header_info)
header_names = tag_names(header_info)

tags = tag_names(headers)
n = n_tags(headers)
; Make lookup table to convert between indices in the header_info structure
; and the headers structure passed into this routine.
tag_id = replicate(99,nnames)
infotag_id = replicate(99,nnames)
for i=0, nnames-1 do begin
  str = header_info.(i)+'_present'
  for j=0,n-1 do if strcmp(str,tags[j],/fold_case) then tag_id[i] = j
  for j=0,n-1 do if strcmp(header_names[i],tags[j],/fold_case) then infotag_id[i] = j
endfor

;put the DWEL header information
for i=0, nnames-1 do begin
  if (tag_id[i] ne 99 and infotag_id[i] ne 99) then begin
    if (headers.(tag_id[i])) then begin
;   help,headers.(infotag_id[i])
        envi_assign_header_value, fid=in_fid, keyword=header_names[i], $
          value=headers.(infotag_id[i])
    endif
  endif
endfor
envi_write_file_header, in_fid
header_info=0b
infotag_id=0b
tag_id=0b
header_names=0b

return,1b

end

