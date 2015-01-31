function list_dwel_headers,inanc,outanc,err
  compile_opt idl2
  ;
  sfile=30
  err=0
  istat=0
  
  if (strlen(outanc) le 0) then begin
    print,'the output file name is empty!!'
    istat=1
    err=1
    goto,out
  endif
  ;make sure path to output exists (can be used to set it up)
  o_dir = file_dirname(outanc)
  file_mkdir, o_dir
  ;
  ;now open the output file and record some stuff
  ;get the output file set up
  ;first see if a file like this exists and is open in envi
  if (file_test(outanc)) then begin
    fids=envi_get_file_ids()
    if(fids[0] ne -1) then begin
      for i=0,n_elements(fids)-1 do begin
        envi_file_query,fids[i],fname=tname
        if (strtrim(strlowcase(outanc),2) eq $
          strtrim(strlowcase(tname),2)) then begin
          envi_file_mng,id=fids[i],/remove
        endif
      endfor
    endif
  endif
  ;open the output file
  openw,sfile,outanc,/get_lun,error=error
  if (error ne 0) then begin
    print,'error opening outanc!!'
    err=2
    istat=1
    goto,out
  endif
  printf,sfile,strtrim('Running Header List for DWEL info',2)
  printf,sfile,strtrim('Run made at: '+systime(),2)
  printf,sfile,strtrim('Target File='+strtrim(inanc,2),2)
  flush,sfile
  
  text_err=0
  envi_open_file,inanc,r_fid=anc_fid,/no_interactive_query,/no_realize
  if (anc_fid eq -1) then begin
    print,'Processing stopped! Error opening input data file: '+strtrim(inanc,2)
    err=3
    istat=1
    goto,out
  endif
  
  ;get the input image dimensions and other info
  envi_file_query, anc_fid, ns=Nshots, nl=Nscans, nb=nb_anc
  
  ;now get the DWEL headers that are present for the input dwel file
  ;set up a base structure for the  headers
  DWEL_anc_headers={ $
    f_base:inanc $
    }
    
  ;find all of the DWEL headers in the hdr file as defined by FID
  status=DWEL_get_headers(anc_fid,DWEL_anc_headers)
  
  if (~status) then begin
    print,'Processing stopped! Bad FID in DWEL_get_headers for input File'
    envi_file_mng,id=anc_fid,/remove
    err=4
    istat=1
    goto, out
  endif
  
  if ((DWEL_anc_headers.headers_present le 0s) or (~DWEL_anc_headers.run_present)) then begin
    print,'DWEL_anc_headers.headers_present='+strtrim(string(DWEL_anc_headers.headers_present),2)
    print,'DWEL_anc_headers.run_present='+strtrim(string(DWEL_anc_headers.run_present),2)
    print,'Processing stopped! File NOT a valid DWEL Cube file'
    envi_file_mng,id=anc_fid,/remove
    err=5
    istat=1
    goto, out
  endif
  
  title=''
  ;find the title
  match = -1
  for i=0,n_elements(DWEL_anc_headers.DWEL_scan_info)-1 do begin
    if (strmatch(DWEL_anc_headers.DWEL_scan_info[i],'*Scan Description*')) then match=i
  endfor
  if (match ge 0) then begin
    text=strtrim(DWEL_anc_headers.DWEL_scan_info[match],2)
    kk=strpos(text,'=')
    title=strtrim(strmid(text,kk+1),2)
  endif
  printf,sfile,'[Basic Info]'
  printf,sfile,'Scan Description='+strtrim(title,2)
  printf,sfile,'Input File='+strtrim(inanc,2)
  flush,sfile
  
  print,'[Basic Info]'
  print,'Scan Description='+strtrim(title,2)
  print,'Input File='+strtrim(inanc,2)
  
  tagnames=tag_names(dwel_anc_headers)
  print,'number of tags=',n_elements(tagnames)
  
  upper=((n_elements(tagnames)-2)/2)
  print,'upper=',upper
  ku=0
  
  for j=0,upper-1 do begin
    ku=(upper-1)-j
    if (byte(dwel_anc_headers.(2*ku+1))) then begin
      printf,sfile,''
      printf,sfile,'Info Tag ['+strtrim(string(j+1),2)+']='+strtrim(tagnames[2*ku],2)
      ;    print,'Info Tag='+strtrim(tagnames[2*ku],2)+' is Present:'
      for kk=0,n_elements(dwel_anc_headers.(2*ku))-1 do begin
        printf,sfile,strtrim((dwel_anc_headers.(2*ku))[kk],2)
      ;      print,strtrim((dwel_anc_headers.(2*ku))[kk],2)
      endfor
    endif
  endfor
  flush,sfile
  free_lun,sfile,/force
  envi_file_mng,id=anc_fid,/remove
  
  out:
  
  free_lun,sfile,/force
  return,istat
end
