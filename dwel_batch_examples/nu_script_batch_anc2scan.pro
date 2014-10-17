;a script to check encoders and conversions using an ancillary file
pro Nu_Script_Batch_Anc2Scan
  compile_opt idl2
  
  tfile=30
  batch_mode=0b
  bale_out=0b
  err_flag=1b
  anc_fid=0
  
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
  print,'Starting Nu_Script_Batch_Anc2Scan'
  print,' '
  
  envi, /restore_base_save_files
  
  help,/memory
  
  ;clean up any fids which are no longer where they were!
  ;ENVI issue that is annoying and leads to confusion
  clean_envi_file_fids
  
  ;Settings section
  ;Set bale_out=1 if you wish to just do the first runs (two bands of one hdf file) to be sure!
  ;set batch_mode=1 if you wish it to run in batch and close IDL at the end
  ;err_flag=1b indicates an error has occured
  ;always set logfile_set to zero here at this time
  bale_out=0b
  batch_mode=0b
  err_flag=0b
  logfile_set=0b
  
  ;====================================================================
  ;set up logfile
  logfile='Y:\DWEL\Data\Processed_data\TumEE05_07082014\Test_zen_tweak\anc2scan_test_angles_TBR_EE_cube_median_logfile.log'
  if (strlen(logfile) le 0) then begin
    print,'the logfile name is empty!!'
    istat=1
    err_flag=1b
    goto,out
  endif
  ;make sure path to logfile exists (this section can be used to set up path)
  o_dir = file_dirname(logfile)
  file_mkdir, o_dir
  ;inanc is the list of HDF files input
  inanc=[$
    'Y:\DWEL\Data\Processed_data\TumEE05_07082014\TumEE05_waveform_2014-08-07-13-37_1548_cube_ancillary.img' $
    ]
  zen_tweak=[$
    -400L $
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
      print,'input ancillary file does NOT exist!'
      print,'Case='+strtrim(string(k),2)
      print,'File='+strtrim(inanc[k],2)
      err_flag=1b
      goto,out
    endif
    if (strpos(strlowcase(inanc[k]),'ancillary') lt 0) then begin
      print,'input file is NOT an ancillary file!'
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
  printf,tfile,strtrim('Running anc2scan for DWEL info',2)
  printf,tfile,strtrim('Run made at: '+systime(),2)
  printf,tfile,strtrim('Number of cases='+strtrim(string(ncase),2),2)
  flush,tfile
  
  logfile_set=1b
  
  ;write out input ancillary files
  printf,tfile,'Input Ancillary  files:'
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
    printf,tfile,'Input Ancillary file='+strtrim(inanc[k],2)
    
    print,'Outer Run Case='+strtrim(string(k+1),2)
    
    inbase=file_basename(inanc[k])
    nbuf=strpos(inbase,'.img')
    baseout=strtrim(strmid(inbase,0,nbuf),2)
    print,'baseout='+baseout
    
    text_err=0
    envi_open_file, inanc[k],r_fid=anc_fid,/no_interactive_query,/no_realize
    if (anc_fid eq -1) then begin
      print,'Processing stopped! Error opening ancillary data file '+strtrim(inanc[k],2)
      goto, out
    endif
    
    ;get the input image dimensions and other info
    envi_file_query, anc_fid, ns=Nshots, nl=Nscans, nb=nb_anc
    
    ;now get the DWEL headers that are present for the ancillary file
    ;set up a base structure for the  headers
    DWEL_anc_headers={ $
      f_base:inanc[k] $
      }
      
    ;find all of the DWEL headers in the hdr file as defined by FID
    status=DWEL_get_headers(anc_fid,DWEL_anc_headers)
    
    if (not status) then begin
      print,'Processing stopped! Bad FID in DWEL_get_headers for Ancillary File'
      envi_file_mng,id=anc_fid,/remove
      goto, out
    endif
    
    if (DWEL_anc_headers.headers_present le 0s or not DWEL_anc_headers.run_present) then begin
      print,'DWEL_anc_headers.headers_present='+strtrim(string(DWEL_anc_headers.headers_present),2)
      print,'DWEL_anc_headers.run_present='+strtrim(string(DWEL_anc_headers.run_present),2)
      print,'Processing stopped! Ancillary file NOT a valid DWEL Cube file'
      envi_file_mng,id=anc_fid,/remove
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
    printf,tfile,'[Encoder Info]'
    printf,tfile,'Scan Description='+strtrim(title,2)
    printf,tfile,'Input File='+strtrim(inanc[k],2)
    printf,tfile,'Nshots='+strtrim(string(Nshots),2)
    printf,tfile,'Nscans='+strtrim(string(Nscans),2)
    flush,tfile
    
    print,'processing '+strtrim(inanc[k],2)
    
    ;get the mask
    Mask_all=bytarr(Nshots,Nscans)+1b
    dims=[-1,0,Nshots-1,0,Nscans-1]
    Mask_all=byte(envi_get_data(fid=anc_fid,dims=dims,pos=6))
    help,mask_all
    
    ;get encoder info
    scanenc=fltarr(Nshots,Nscans)
    rotaryenc=fltarr(Nshots,Nscans)
    scanenc=float(envi_get_data(fid=anc_fid,dims=dims,pos=2))
    rotaryenc=float(envi_get_data(fid=anc_fid,dims=dims,pos=3))
    
    ;get the zenith and azimuth
    ;set up a structure and push it onto the heap
    shotzen=fltarr(Nshots,Nscans)
    shotazim=fltarr(Nshots,Nscans)
    shotzen=scanenc
    shotazim=rotaryenc
    sav={ $
      Nshots:Nshots,$
      Nscans:Nscans,$
      ShotZen:ShotZen,$
      ShotAzim:ShotAzim $
      }
    ;now put the data on the heap with a pointer
    p_stat=ptr_new(sav,/no_copy)
    status = DWEL_set_theta_phi_oz(p_stat,zen_tweak[k])
    ;put the results into the local arrays
    ShotZen=(*p_stat).ShotZen
    ShotAzim=(*p_stat).ShotAzim
    ptr_free, p_stat
    
    ;now go and do it
    print,'nshots='+strtrim(string(Nshots),2)
    print,'nscans='+strtrim(string(Nscans),2)
    
    ;first the encoders
    min_scan_vec=fltarr(Nshots)
    max_scan_vec=fltarr(Nshots)
    mean_scan_vec=fltarr(Nshots)
    
    for j=0L,Nshots-1L do begin
      pos=where(reform(mask_all[j,*]) ne 0,npos)
      if (npos gt 0) then begin
        min_scan_vec[j]=min(reform(scanenc[j,pos]))
        max_scan_vec[j]=max(reform(scanenc[j,pos]))
        ;    mean_scan_vec[j]=total(reform(scanenc[j,pos]))/float(npos)
        mean_scan_vec[j]=median(reform(scanenc[j,pos]))
      endif else begin
        min_scan_vec[j]=0.0
        max_scan_vec[j]=0.0
        mean_scan_vec[j]=0.0
      endelse
      pos=0b
    endfor
    
    min_rot_vec=fltarr(Nscans)
    max_rot_vec=fltarr(Nscans)
    mean_rot_vec=fltarr(Nscans)
    
    for j=0L,Nscans-1L do begin
      pos=where(reform(mask_all[*,j]) ne 0,npos)
      if (npos gt 0) then begin
        min_rot_vec[j]=min(reform(rotaryenc[pos,j]))
        max_rot_vec[j]=max(reform(rotaryenc[pos,j]))
        ;    mean_rot_vec[j]=total(reform(rotaryenc[pos,j]))/float(npos)
        mean_rot_vec[j]=median(reform(rotaryenc[pos,j]))
      endif else begin
        min_rot_vec[j]=0.0
        max_rot_vec[j]=0.0
        mean_rot_vec[j]=0.0
      endelse
      pos=0b
    endfor
    
    ;now the angles
    min_zen_vec=fltarr(Nshots)
    max_zen_vec=fltarr(Nshots)
    mean_zen_vec=fltarr(Nshots)
    
    for j=0L,Nshots-1L do begin
      pos=where(reform(mask_all[j,*]) ne 0,npos)
      if (npos gt 0) then begin
        min_zen_vec[j]=min(reform(shotzen[j,pos]))
        max_zen_vec[j]=max(reform(shotzen[j,pos]))
        ;    mean_zen_vec[j]=total(reform(shotzen[j,pos]))/float(npos)
        mean_zen_vec[j]=median(reform(shotzen[j,pos]))
      endif else begin
        min_zen_vec[j]=0.0
        max_zen_vec[j]=0.0
        mean_zen_vec[j]=0.0
      endelse
      pos=0b
    endfor
    
    min_azim_vec=fltarr(Nscans)
    max_azim_vec=fltarr(Nscans)
    mean_azim_vec=fltarr(Nscans)
    
    for j=0L,Nscans-1L do begin
      pos=where(reform(mask_all[*,j]) ne 0,npos)
      if (npos gt 0) then begin
        min_azim_vec[j]=min(reform(shotazim[pos,j]))
        max_azim_vec[j]=max(reform(shotazim[pos,j]))
        ;    mean_azim_vec[j]=total(reform(shotazim[pos,j]))/float(npos)
        mean_azim_vec[j]=median(reform(shotazim[pos,j]))
      endif else begin
        min_azim_vec[j]=0.0
        max_azim_vec[j]=0.0
        mean_azim_vec[j]=0.0
      endelse
      pos=0b
    endfor
    
    buf=''
    
    printf,tfile,'Scan Zenith Encoder Info, i,MinScan,MaxScan,MedScan,MinZen,MaxZen,MedZen'
    for j=0L,Nshots-1L do begin
      buf=strtrim(string(j),2)+','+strtrim(string(min_scan_vec[j]),2) $
        +','+strtrim(string(max_scan_vec[j]),2) $
        +','+strtrim(string(mean_scan_vec[j]),2) $
        +','+strtrim(string(min_zen_vec[j]),2) $
        +','+strtrim(string(max_zen_vec[j]),2) $
        +','+strtrim(string(mean_zen_vec[j]),2)
      printf,tfile,strtrim(buf,2)
      buf=''
    endfor
    
    printf,tfile,'Rotary Azimuth Encoder Info, i,MinRot,MaxRot,MedRot,MinAzim,MaxAzim,MedAzim'
    for j=0L,Nscans-1L do begin
      buf=strtrim(string(j),2)+','+strtrim(string(min_rot_vec[j]),2) $
        +','+strtrim(string(max_rot_vec[j]),2) $
        +','+strtrim(string(mean_rot_vec[j]),2) $
        +','+strtrim(string(min_azim_vec[j]),2) $
        +','+strtrim(string(max_azim_vec[j]),2) $
        +','+strtrim(string(mean_azim_vec[j]),2)
      printf,tfile,strtrim(buf,2)
      buf=''
    endfor
    
    rotaryenc=0b
    scanenc=0b
    mask=0b
    envi_file_mng,id=anc_fid,/remove
    ;check open files
    fids=envi_get_file_ids()
    if(fids[0] ge 0) then begin
      for i=0,n_elements(fids)-1 do begin
        envi_file_query,fids[i],fname=tname
        print,'Open File '+strtrim(string(i+1),2)+'='+strtrim(tname,2)
        envi_file_mng,id=fids[i],/remove
      endfor
    endif
    
    print,'End of outer run case='+strtrim(string(ncase+1),2)
    help,/memory
    
  endfor
  
  T_total=systime(1)
  print,'TOTAL Elapsed Time for all Cases was '+strtrim(string(float(T_total-T_outer)/60.0,format='(f12.3)'),2)+' minutes'
  
  printf,tfile,'TOTAL Elapsed Time for all Cases was '+strtrim(string(float(T_total-T_outer)/60.0,format='(f12.3)'),2)+' minutes'
  printf,tfile,'Nu_Script_Batch_Anc2Scan finished'
  flush,tfile
  free_lun,tfile,/force
  logfile_set=0b
  
  out:
  free_lun,tfile,/force
  logfile_set=0b
  envi_file_mng,id=anc_fid,/remove
  ;
  print,'Nu_Script_Batch_Anc2Scan finished'
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
