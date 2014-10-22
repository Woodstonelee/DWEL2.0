function dwel_get_config_info,config_file,Config_Info,consum
  compile_opt idl2
  
  inlun=30
  err=0
  consum=''
  ;
  ;get configuration info
  if (~file_test(Config_File)) then return,1
  openr,inlun,config_file,/get_lun,error=err
  if (err ne 0) then begin
    print,'error opening configuration file'
    print,'error number='+strtrim(string(err),2)
    print,'File='+strtrim(config_file,2)
    return,1
  endif
  buf=''
  nrec=0
  reading:
  if(eof(inlun)) then goto,done
  readf,inlun,buf
  buf=strtrim(buf,2)
  ;
  if(buf eq '') then goto,reading
  npl=strpos(buf,' = ')
  if (npl gt -1) then begin
    left=strmid(buf,0,npl)
    right=strmid(buf,npl+3)
    buf=strtrim(left,2)+'='+strtrim(right,2)
  endif
  npl=strpos(buf,'comment=')
  if (npl gt -1) then begin
    consum=strtrim(strcompress(strmid(buf,npl+8)),2)
    if (nrec le 0) then con_entries=[buf] else con_entries=[con_entries,buf]
    nrec=nrec+1
  endif 
  goto,reading
  done:
  
  free_lun,inlun,/force
  
  for k=0,nrec-1 do begin
    con_entries[k]=strtrim(con_entries[k],2)
  endfor
  Config_Info=[con_entries]
  ;
  return,0
end
