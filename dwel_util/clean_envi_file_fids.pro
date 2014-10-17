pro clean_envi_file_fids
  compile_opt idl2
  ;
  ;routine to clean up files with fids not cannot be found
  ;with their current stored path & name
  ;
  fids=envi_get_file_ids()
  if(fids[0] ne -1) then begin
    for i=0,n_elements(fids)-1 do begin
      envi_file_query,fids[i],fname=tname
      if (~file_test(tname)) then begin
        envi_file_mng,id=fids[i],/remove
      endif
    endfor
  endif
  
  return
  
end
