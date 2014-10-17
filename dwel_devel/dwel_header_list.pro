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
  
  if (not status) then begin
    print,'Processing stopped! Bad FID in DWEL_get_headers for input File'
    envi_file_mng,id=anc_fid,/remove
    err=4
    istat=1
    goto, out
  endif
  
  if (DWEL_anc_headers.headers_present le 0s or not DWEL_anc_headers.run_present) then begin
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

;a script to List DWEL headers from files
pro DWEL_header_List
  compile_opt idl2
  
  tfile=30
  batch_mode=0b
  bale_out=0b
  err_flag=1b
  anc_fid=0
  logfile_set=0b
  logfile=''
  
  day=['Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday']
  month=['January','February','March','April','May','June','July', $
    'August','September','October','November','December']
    
  ;Now setup protective aunty-Catch to watch out for errors
  error_status=0
  catch, error_status
  if (error_status ne 0) then begin
    catch,/cancel
    help,/last,output=out
    info_text=[strtrim(string('Catch trapped an Error !'),2),$
      strtrim('Error Name : '+strtrim(!err_string,2),2),$
      strtrim(string('Error Number: ',error_status,$
      format='(a,i5)'),2),$
      strtrim('Last Message: '+strtrim(out,2),2)]
    ;ie write this out to the IDL console when possible
    for j=0,n_elements(info_text)-1 do begin
      print,info_text[j]
    endfor
    if (logfile_set) then begin
      for j=0,n_elements(info_text)-1 do begin
        printf,tfile,info_text[j]
      endfor
      flush,tfile
    endif
    err_flag=1b
    goto, out
  endif
  
  print,' '
  print,'Starting DWEL_header_List'
  print,' '
  
  envi, /restore_base_save_files
  
  help,/memory
  
  ;clean up any fids which are no longer where they were!
  ;ENVI issue that is annoying and leads to confusion
  clean_envi_file_fids
  
  ;====================================================================
  ;Settings section
  ;Set bale_out=1 if you wish to just do the first runs (two bands of one hdf file) to be sure!
  ;set batch=1 if you wish it to run in batch and close IDL at the end
  ;err_flag=1b indicates an error has occured
  ;always set logfile_set to zero here at this time
  bale_out=0b
  batch_mode=0b
  err_flag=0b
  logfile_set=0b
  ;set up logfile
  logfile='Y:\DWEL\Data\Processed_data\TumEE05_07082014\dwel_header_list_test_logfile.log'
  if (strlen(logfile) le 0) then begin
    print,'the logfile name is empty!!'
    istat=1
    err_flag=1b
    goto,out
  endif
  ;make sure path to logfile exists (can be used to set it up)
  o_dir = file_dirname(logfile)
  file_mkdir, o_dir
  ;inanc is the list of DWEL files input
  inanc=[$
    ;'Y:\DWEL\Data\Processed_data\TumEE05_07082014\TumEE05_waveform_2014-08-07-13-37_1064_cube_pulse_filt_update_at_proj_full.img',$
    'Y:\DWEL\Data\Processed_data\TumEE05_07082014\TumEE05_waveform_2014-08-07-13-37_1064_cube_pulse_filt_update_atp_ptcl_pcinfo.img' $
    ]
    
  ;====================================================================
  ;now start to get things checked
  ncase=n_elements(inanc)
  
  ;print,'Number of cases=',ncase
  if (ncase le 0) then begin
    print,'Number of cases invalid!'
    err_flag=1b
    goto,out
  endif
  
  print,'Number of Cases='+strtrim(string(ncase),2)
  
  for k=0,ncase-1 do begin
    inanc[k]=strtrim(inanc[k],2)
    if (~file_test(inanc[k])) then begin
      print,'input DWEL file does NOT exist!'
      print,'Case='+strtrim(string(k),2)
      print,'File='+strtrim(inanc[k],2)
      err_flag=1b
      goto,out
    endif
  endfor
  
  print,'input checking complete'
  
  ;now open the log file and record some stuff
  ;get the log file set up
  ;first see if a file like this exists and is open in envi
  if (file_test(logfile)) then begin
    fids=envi_get_file_ids()
    if(fids[0] ne -1) then begin
      for i=0,n_elements(fids)-1 do begin
        envi_file_query,fids[i],fname=tname
        if (strtrim(strlowcase(logfile),2) eq $
          strtrim(strlowcase(tname),2)) then begin
          envi_file_mng,id=fids[i],/remove
        endif
      endfor
    endif
  endif
  ;open the log file
  openw,tfile,logfile,/get_lun,error=error
  if (error ne 0) then begin
    print,'error opening logfile!!'
    err_flag=1b
    goto,out
  endif
  printf,tfile,strtrim('Running Header List for DWEL info',2)
  printf,tfile,strtrim('Run made at: '+systime(),2)
  printf,tfile,strtrim('Number of cases='+strtrim(string(ncase),2),2)
  flush,tfile
  
  logfile_set=1b
  
  ;write out input files
  printf,tfile,'Input files:'
  for j=0L,ncase-1L do begin
    printf,tfile,strtrim(inanc[j],2)
  endfor
  
  ;print,buf
  buf=''
  
  ;all has been set up and checked as far as possible
  printf,tfile,'All Set and Ready to run for '+strtrim(string(ncase),2)+' runs'
  flush,tfile
  
  ;envi_batch_init if batch_mode is set to 1 (true)
  if (batch_mode) then begin
    envi_batch_init,/NO_STATUS_WINDOW
    envi_batch_status_window,/off
  endif
  
  T_outer=systime(1)
  
  ;Main loop over the cases
  for k=0,ncase-1 do begin
  
    printf,tfile,'Outer Run Case='+strtrim(string(k+1),2)
    printf,tfile,'Input file='+strtrim(inanc[k],2)
    
    print,'Outer Run Case='+strtrim(string(k+1),2)
    
    inbase=file_basename(inanc[k])
    nbuf=strpos(inbase,'.img')
    baseout=strtrim(strmid(inbase,0,nbuf),2)
    print,'baseout='+baseout
    
    outanc=strtrim(strmid(inanc[k],0,strpos(inanc[k], '.', /reverse_search))+'_Headers_List.txt',2)
    
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
    
    istat=list_dwel_headers(inanc[k],outanc,err)
    
    if (istat or (err gt 0)) then begin
      print,'return from list_dwel_headers with error'
      print,'error number='+strtrim(string(err),2)
      err_flag=1b
      goto,out
    endif
    
    ;check open files
    fids=envi_get_file_ids()
    if(fids[0] ge 0) then begin
      for i=0,n_elements(fids)-1 do begin
        envi_file_query,fids[i],fname=tname
        print,'Open File '+strtrim(string(i+1),2)+'='+strtrim(tname,2)
        envi_file_mng,id=fids[i],/remove
      endfor
    endif
    
    print,'End of outer run case='+strtrim(string(k+1),2)
    help,/memory
    
  endfor
  
  T_total=systime(1)
  print,'TOTAL Elapsed Time for all Cases was '+strtrim(string(float(T_total-T_outer)/60.0,format='(f12.3)'),2)+' minutes'
  
  printf,tfile,'TOTAL Elapsed Time for all Cases was '+strtrim(string(float(T_total-T_outer)/60.0,format='(f12.3)'),2)+' minutes'
  printf,tfile,'DWEL_header_List finished'
  flush,tfile
  free_lun,tfile,/force
  logfile_set=0b
  
  out:
  free_lun,tfile,/force
  logfile_set=0b
  envi_file_mng,id=anc_fid,/remove
  ;
  print,'DWEL_header_List finished'
  if (err_flag) then begin
    print,'An ERROR has occurred! Check the Log File'
    print,'Log File Name='+strtrim(logfile,2)
  endif
  heap_gc,/verbose
  if (batch_mode and (err_flag eq 0) and (bale_out eq 0)) then begin
    envi_batch_exit,/exit_idl,/no_confirm
  endif
;
end
