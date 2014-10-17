; basefix and satfix and then Filter the basefix and sat fixed data and then run new filtered basefix
pro Nu_Script_Batch_pulse_filt
  compile_opt idl2
  ;
  bale_out=0b
  batch_mode=1b
  err_flag=0b
  logfile_set=0b
  tfile=30
  get_info_stats=1b
  
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
  print,'Starting Nu_Script_Batch_pulse_filt'
  print,' '
  
  envi,/restore_base_save_files
  help,/memory
  
  ;clean up any fids which are no longer where they were!
  ;ENVI issue that is annoying and leads to confusion
  clean_envi_file_fids
  
  ;set up for logfile
  logfile_set=0b
  err_flag=0b
  
  ;between the ==== lines are the settings - can be saved to ascii file for
  ;convenience, re-runs, templates etc etc
  ;====================================================================
  ;Settings for nu_script_batch_pulse_filt
  ;first set the log file name
  logfile='Y:\DWEL\Data\Processed_data\TumEE05_07082014\pulse_filtfix_debug_last2_logfile.log'
  
  ;incfile is the file to be filtered etc
  incfile=[$
    ;'Y:\DWEL\Data\Processed_data\wire10092014\wire10_waveform_2014-09-10-13-44_1064_cube.img',$
    ;'Y:\DWEL\Data\Processed_data\wire10092014\wire10_waveform_2014-09-10-13-44_1548_cube.img',$
    ;'Y:\DWEL\Data\Processed_data\BE004_03092014_Ptcl_Test\BE004_waveform_2014-09-03-13-37_1064_cube.img',$
    ;'Y:\DWEL\Data\Processed_data\BE004_03092014_Ptcl_Test\BE004_waveform_2014-09-03-13-37_1548_cube.img',$
    'Y:\DWEL\Data\Processed_data\TumEE05_07082014\TumEE05_waveform_2014-08-07-13-37_1064_cube.img',$
    'Y:\DWEL\Data\Processed_data\TumEE05_07082014\TumEE05_waveform_2014-08-07-13-37_1548_cube.img',$
    'Y:\DWEL\Data\Processed_data\Pye2_Ptcld\pye2a_waveform_2014-06-21-10-55_1064_cube.img',$
    'Y:\DWEL\Data\Processed_data\Pye2_Ptcld\pye2a_waveform_2014-06-21-10-55_1548_cube.img' $
    ]
  ancfile=[$
    ;'Y:\DWEL\Data\Processed_data\wire10092014\wire10_waveform_2014-09-10-13-44_1064_cube_ancillary.img',$
    ;'Y:\DWEL\Data\Processed_data\wire10092014\wire10_waveform_2014-09-10-13-44_1548_cube_ancillary.img',$
    ;'Y:\DWEL\Data\Processed_data\BE004_03092014_Ptcl_Test\BE004_waveform_2014-09-03-13-37_1064_cube_ancillary.img',$
    ;'Y:\DWEL\Data\Processed_data\BE004_03092014_Ptcl_Test\BE004_waveform_2014-09-03-13-37_1548_cube_ancillary.img',$
    'Y:\DWEL\Data\Processed_data\TumEE05_07082014\TumEE05_waveform_2014-08-07-13-37_1064_cube_ancillary.img',$
    'Y:\DWEL\Data\Processed_data\TumEE05_07082014\TumEE05_waveform_2014-08-07-13-37_1548_cube_ancillary.img',$
    'Y:\DWEL\Data\Processed_data\Pye2_Ptcld\pye2a_waveform_2014-06-21-10-55_1064_cube_ancillary.img',$
    'Y:\DWEL\Data\Processed_data\Pye2_Ptcld\pye2a_waveform_2014-06-21-10-55_1548_cube_ancillary.img' $
    ]
  ;outpath is the path for the output files - names set up automatically
  outpath=[$
    ;'Y:\DWEL\Data\Processed_data\wire10092014\',$
    ;'Y:\DWEL\Data\Processed_data\wire10092014\',$
    ;'Y:\DWEL\Data\Processed_data\BE004_03092014_Ptcl_Test\',$
    ;'Y:\DWEL\Data\Processed_data\BE004_03092014_Ptcl_Test\',$
    'Y:\DWEL\Data\Processed_data\TumEE05_07082014\',$
    'Y:\DWEL\Data\Processed_data\TumEE05_07082014\',$
    'Y:\DWEL\Data\Processed_data\Pye2_Ptcld\',$
    'Y:\DWEL\Data\Processed_data\Pye2_Ptcld\' $
    ]
  zen_tweak=[$
    ;-1200L,-1200L $ ;saturated Clab wire run
    ;,-700L,-700L $ ;This is the Vinyard BE004
    -400L,-400L $ ;This is TumEE05
    ,3154,3154 $ ;This is pye2
    ]
  ;
  ;settings for flow control
  bale_out=0b
  batch_mode=1b
  ;other settings
  ;casing range should be the reflectance panel - used in basefix
  Casing_Range=[170.0,180.0]
  ;settings for anc2at
  max_zenith_angle_at=117.0
  output_resolution_at=2.5
  ;
  ;====================================================================
  ;set up logfile
  if (strlen(logfile) le 0) then begin
    print,'the logfile name is empty!!'
    istat=1
    err_flag=1b
    goto,out
  endif
  ;make sure path to logfile exists (can be used to set it up)
  o_dir = file_dirname(logfile)
  file_mkdir, o_dir
  
  ;now do basic testing of the input
  incfile=strtrim(incfile,2)
  ncase=n_elements(incfile)
  if ((ncase le 0) or (n_elements(outpath) ne ncase) or (n_elements(zen_tweak) ne ncase)) then begin
    print,'Number of cases mismatch in one or more of the input arrays of info!'
    print,'number of cases='+strtrim(string(ncase),2)
    err_flag=1b
    goto,out
  endif
  
  print,'ncase='+strtrim(string(ncase),2)
  
  inwl=intarr(ncase)
  
  for k=0,ncase-1 do begin
    incfile[k]=strtrim(incfile[k],2)
    if (~file_test(incfile[k])) then begin
      print,'an input dwel cube file does NOT exist!'
      print,'Case='+strtrim(string(k),2)
      print,'File='+strtrim(incfile[k],2)
      err_flag=1b
      goto,out
    endif
    ancfile[k]=strtrim(ancfile[k],2)
    if (~file_test(ancfile[k])) then begin
      print,'the ancillary dwel cube file does NOT exist!'
      print,'Case='+strtrim(string(k),2)
      print,'File='+strtrim(ancfile[k],2)
      err_flag=1b
      goto,out
    endif
    ;check each input file is a DWEL file ready for base and sat fix
    ;Open DWEL cube file
    envi_open_file, incfile[k], r_fid=infile_fid,/no_realize,/no_interactive_query
    if (infile_fid eq -1) then begin
      print,strtrim('Error opening input file at case='+strtrim(string(k),2),2)
      print,'File='+strtrim(incfile[k],2)
      err_flag=1b
      goto,out
    endif
    envi_file_query, infile_fid, ns=ns, nl=nl, nb=nb, wl=wl, $
      xstart=xstart, ystart=ystart, data_type=type, $
      interleave=ftype, fname=fname, dims=dims
    f_base=file_basename(fname)
    ;now get the DWEL headers that are present
    ;set up a base structure for the DWEL headers
    DWEL_headers={ $
      f_base:f_base $
      }
    ;find all of the DWEL headers in the hdr file
    status=DWEL_get_headers(infile_fid,DWEL_headers)
    if (not status) then begin
      print,strtrim('Bad call to DWEL_get_headers!!',2)
      print,'Case='+strtrim(string(k),2)
      print,'File='+strtrim(incfile[k],2)
      envi_file_mng,id=infile_fid,/remove
      err_flag=1b
      goto,out
    endif
    if (DWEL_headers.headers_present le 0s or not DWEL_headers.run_present) then begin
      print,strtrim('Input file is NOT a DWEL file',2)
      print,'Case='+strtrim(string(k),2)
      print,'File='+strtrim(incfile[k],2)
      envi_file_mng,id=infile_fid,/remove
      err_flag=1b
      goto,out
    endif
    info=DWEL_headers.dwel_adaptation
    ;now get the DWEL wavelength
    match = -1
    for i=0,n_elements(info)-1 do begin
      if (strmatch(info[i],'*Wavelength=*', /fold_case)) then match=i
    endfor
    if match ge 0 then begin
      text=strtrim(info[match],2)
      kpos=strpos(text,'=')
      wavelength=fix(strtrim(strmid(text,kpos+1,4),2))
    endif else begin
      print,strtrim('Input file does NOT have a DWEL file wavelength',2)
      print,'Case='+strtrim(string(k),2)
      print,'File='+strtrim(incfile[k],2)
      envi_file_mng,id=infile_fid,/remove
      err_flag=1b
      goto,out
    endelse
    inwl[k]=wavelength
    envi_file_mng,id=infile_fid,/remove
    wl_string=strtrim(string(wavelength),2)
    if (strpos(incfile[k],wl_string) lt 0) then begin
      print,'Warning - Wavelength from file does not appear in name'
      print,'Wavelength from Header='+wl_string
      print,'File='+strtrim(incfile[k],2)
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
  printf,tfile,strtrim('Running base_satfix, filter, refresh basefix plus anc2at for DWEL pre-processing',2)
  printf,tfile,strtrim('Run made at: '+systime(),2)
  printf,tfile,strtrim('Number of cases='+strtrim(string(ncase),2),2)
  flush,tfile
  
  logfile_set=1b
  
  ;write out input files
  printf,tfile,'Input files:'
  for j=0L,ncase-1L do begin
    printf,tfile,strtrim(incfile[j],2)
  endfor
  ;write out input ancillary files
  printf,tfile,'Input ancillary files:'
  for j=0L,ncase-1L do begin
    printf,tfile,strtrim(ancfile[j],2)
  endfor
  ;write output file paths
  printf,tfile,'Output file paths:'
  for j=0L,ncase-1L do begin
    printf,tfile,strtrim(outpath[j],2)
  endfor
  ;now list the wavelengths
  buf=strtrim(string(inwl),2)
  buf='['+strtrim(strjoin(buf,','),2)+']'
  printf,tfile,'Input wavelengths='+buf
  ;now list the zen_tweaks
  buf=strtrim(string(zen_tweak),2)
  buf='['+strtrim(strjoin(buf,','),2)+']'
  printf,tfile,'Input Zenith Tweaks='+buf
  
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
  
  for k=0,ncase-1 do begin
  
    print,''
    printf,tfile,'Outer Run Case='+strtrim(string(k+1),2)
    printf,tfile,'Input Cube file='+strtrim(incfile[k],2)
    printf,tfile,'Input Ancillary file='+strtrim(ancfile[k],2)
    
    inbase=file_basename(incfile[k])
    nbuf=strpos(inbase,'.img')
    baseout=strtrim(strmid(inbase,0,nbuf),2)
    print,'baseout='+baseout
    
    ;
    outdir=outpath[k]
    incube=strtrim(incfile[k],2)
    
    BasefixImageFile=outdir+baseout+'_bsfix_image.img'
    outcube_filter=outdir+baseout+'_pulse_filt.img'
    outupdatedfile=outdir+baseout+'_pulse_filt_update.img'
    DWEL_AT_File=outdir+baseout+'_cube_Anc2AT_project.img'
    
    print,''
    print,'outcube_filter='+outcube_filter
    
    ; Create the output directory (if it already exists, this code will do nothing)
    o_dir = file_dirname(BasefixImageFile)
    file_mkdir, o_dir
    ; Create the output directory (if it already exists, this code will do nothing)
    o_dir = file_dirname(outcube_filter)
    file_mkdir, o_dir
    ; Create the output directory (if it already exists, this code will do nothing)
    o_dir = file_dirname(outupdatedfile)
    file_mkdir, o_dir
    ; Create the output directory (if it already exists, this code will do nothing)
    o_dir = file_dirname(DWEL_AT_File)
    file_mkdir, o_dir
    
    ;check files are not open in envi
    if(file_test(BasefixImageFile)) then begin
      fids=envi_get_file_ids()
      if(fids[0] eq -1) then begin
      ;
      endif else begin
        for i=0,n_elements(fids)-1 do begin
          envi_file_query,fids[i],fname=tname
          if (strtrim(strlowcase(BasefixImageFile),2) eq $
            strtrim(strlowcase(tname),2)) then begin
            envi_file_mng,id=fids[i],/remove
            print,'BasefixImageFile fid removed'
          endif
        endfor
      endelse
    endif
    ;
    if(file_test(outcube_filter)) then begin
      fids=envi_get_file_ids()
      if(fids[0] eq -1) then begin
      ;
      endif else begin
        for i=0,n_elements(fids)-1 do begin
          envi_file_query,fids[i],fname=tname
          if (strtrim(strlowcase(outcube_filter),2) eq $
            strtrim(strlowcase(tname),2)) then begin
            envi_file_mng,id=fids[i],/remove
            print,'outcube_filter fid removed'
          endif
        endfor
      endelse
    endif
    ;
    if(file_test(outupdatedfile)) then begin
      fids=envi_get_file_ids()
      if(fids[0] eq -1) then begin
      ;
      endif else begin
        for i=0,n_elements(fids)-1 do begin
          envi_file_query,fids[i],fname=tname
          if (strtrim(strlowcase(outupdatedfile),2) eq $
            strtrim(strlowcase(tname),2)) then begin
            envi_file_mng,id=fids[i],/remove
            print,'outupdatedfile fid removed'
          endif
        endfor
      endelse
    endif
    ;
    if(file_test(DWEL_AT_File)) then begin
      fids=envi_get_file_ids()
      if(fids[0] eq -1) then begin
      ;
      endif else begin
        for i=0,n_elements(fids)-1 do begin
          envi_file_query,fids[i],fname=tname
          if (strtrim(strlowcase(DWEL_AT_File),2) eq $
            strtrim(strlowcase(tname),2)) then begin
            envi_file_mng,id=fids[i],/remove
            print,'DWEL_AT_File fid removed'
          endif
        endfor
      endelse
    endif
    ;
    
    T_start=systime(1)
    
    ; next run gets baseline and satfixed data
    ;
    
    ; Create the output directory (if it already exists, this code will do nothing)
    o_dir = file_dirname(BasefixImageFile)
    file_mkdir, o_dir
    
    print,'Calling DWEL_baseline_sat_fix_cmd_oz'
    
    buf=''
    buf=strtrim(string(Casing_Range),2)
    buf='['+strtrim(strjoin(buf,','),2)+']'
    printf,tfile,'Casing Range(deg) ='+buf
    
    err=0
    DWEL_baseline_sat_fix_cmd_oz, incfile[k],Ancfile[k],BasefixImageFile,Casing_Range,get_info_stats,zen_tweak[k],err
    
    if (err ne 0) then begin
      printf,tfile,'Error in call to DWEL_baseline_sat_fix_cmd_oz from Nu_Script_Batch_pulse_filt'
      printf,tfile,'Local error code='+strtrim(string(err),2)
      err_flag=1b
      goto,out
    endif
    
    if(~file_test(BasefixImageFile)) then begin
      printf,tfile,'Output file from DWEL_baseline_sat_fix_cmd_oz does NOT exist!'
      printf,tfile,'expected file is='+strtrim(BasefixImageFile,2)
      goto,out
      err_flag=1b
    endif
    ;
    ;Set up expected ancillary file
    
    AncillaryFile = strtrim(strmid(BasefixImageFile,0,strpos(BasefixImageFile, '.', /reverse_search))+'_ancillary.img',2)
    if(~file_test(AncillaryFile)) then begin
      printf,tfile,'Expected Output ancillary file from DWEL_baseline_sat_fix_cmd_oz does NOT exist!'
      printf,tfile,'expected file is='+strtrim(AncillaryFile,2)
      goto,out
      err_flag=1b
    endif
    ;
    
    print,'Output Basefix Image='+strtrim(BasefixImageFile,2)
    
    T_first=systime(1)
    
    ;swop the filters so that data are filtered with other wavelength pulse
    wavelength=inwl[k]
    sel_wl=wavelength
    if (wavelength eq 1064) then sel_wl=1548 else sel_wl=1064
    
    DWEL_pulse_model_dual_oz,sel_wl,i_val,t_val,r_val,p_range,p_time,pulse,t_fwhm,r_fwhm
    
    print,''
    print,'Number of values in filter='+strtrim(string(n_elements(pulse)),2)
    
    pulse=pulse/total(pulse)
    ierr=0
    maxval=max(pulse,mpos)
    
    print,'Mpos='+strtrim(string(mpos),2)
    ;run a general filter (not dwel specific)
    istat=dwel_general_filter(BasefixImageFile,pulse,mpos,outcube_filter,ierr)
    
    if(~file_test(outcube_filter) or (ierr ne 0) or (istat ne 0)) then begin
      print,'Output file from DWEL_general_filter does NOT exist or error!'
      print,'expected file is='+strtrim(outcube_filter,2)
      if (ierr ne 0) then print,'Local error code='+strtrim(string(ierr),2)
      printf,tfile,'Output file from DWEL_general_filter does NOT exist or error!'
      printf,tfile,'expected file is='+strtrim(outcube_filter,2)
      if (ierr ne 0) then printf,tfile,'Local error code='+strtrim(string(ierr),2)
      err_flag=1b
      goto,out
    endif
    
    print,'Output Filtered Image='+strtrim(outcube_filter,2)
    printf,tfile,'Output Filtered Image='+strtrim(outcube_filter,2)
    
    T_sec=systime(1)
    
    ;now update all the stats and ancillary file etc
    
    ;new routine to finalise all the calibrations, Tzero etc
    ;note ancillary file is from previous run but file is filtered by the general filter
    DWEL_Filtered_FixBase_Cmd_oz, outcube_filter, AncillaryFile, outupdatedfile, get_info_stats, zen_tweak[k], err
    
    if (err gt 0) then begin
      print,'DWEL_Filtered_FixBase_Cmd returned with error'
      print,'Local error number='+strtrim(string(err),2)
      printf,tfile,'DWEL_Filtered_FixBase_Cmd returned with error'
      printf,tfile,'Local error number='+strtrim(string(err),2)
      err_flag=1b
      goto,out
    endif
    
    if(~file_test(outupdatedfile)) then begin
      print,'Output updated file does NOT exist or error!'
      print,'expected file is='+strtrim(outupdatedfile,2)
      printf,tfile,'Output updated file does NOT exist or error!'
      printf,tfile,'expected file is='+strtrim(outupdatedfile,2)
      err_flag=1b
      goto,out
    endif
    
    print,'Output updated Filtered Image='+strtrim(outupdatedfile,2)
    printf,tfile,'Output updated Filtered Image='+strtrim(outupdatedfile,2)
    
    ;Set up expected ancillary file
    
    AncillaryFile = strtrim(strmid(outupdatedfile,0,strpos(outupdatedfile,'.', /reverse_search))+'_ancillary.img',2)
    if(~file_test(AncillaryFile)) then begin
      print,'Expected Output ancillary file from DWEL_baseline_sat_fix_cmd_oz does NOT exist!'
      print,'expected file is='+strtrim(AncillaryFile,2)
      goto,out
      err_flag=1b
    endif
    ;
    
    print,'Output ancillary updated Filtered Image file created: '+strtrim(AncillaryFile,2)
    printf,tfile,'Output ancillary updated Filtered Image file created: '+strtrim(AncillaryFile,2)
    
    ;
    ; next run is to get at projected qlook
    
    T_third=systime(1)
    
    printf,tfile,'Max Zenith Angle='+strtrim(string(max_zenith_angle_at),2)+' deg'
    printf,tfile,'Output Resolution='+strtrim(string(output_resolution_at),2)+' mrad'
    
    print,'Calling dwel_anc2at'
    
    err=0
    ;get a Qlook at projected image to see that all is well
    dwel_anc2at_oz, AncillaryFile, DWEL_AT_File, max_zenith_angle_at, output_resolution_at,zen_tweak[k],err
    
    if (err gt 0) then begin
      print,'dwel_anc2at returned with error'
      print,'Local error number='+strtrim(string(err),2)
      printf,tfile,'dwel_anc2at returned with error'
      printf,tfile,'Local error number='+strtrim(string(err),2)
      err_flag=1b
      goto,out
    endif
    
    if(~file_test(DWEL_AT_File)) then begin
      printf,tfile,'Output file from dwel_anc2at does NOT exist!'
      printf,tfile,'expected file is='+strtrim(DWEL_AT_File,2)
      err_flag=1b
      goto,out
    endif
    
    printf,tfile,'Output AT Qlook file='+strtrim(DWEL_AT_File,2)
    
    T_end=systime(1)
    print,' '
    print,'Time for Base & Sat Fix was '+strtrim(string(float(T_first-T_start)/60.0,format='(f12.3)'),2)+' minutes'
    print,'Time for Filter was '+strtrim(string(float(T_sec-T_first)/60.0,format='(f12.3)'),2)+' minutes'
    print,'Time for Base Re-Fix was '+strtrim(string(float(T_third-T_sec)/60.0,format='(f12.3)'),2)+' minutes'
    print,'Time for anc2at was '+strtrim(string(float(T_end-T_third)/60.0,format='(f12.3)'),2)+' minutes'
    print,'Total Elapsed Time was '+strtrim(string(float(T_end-T_start)/60.0,format='(f12.3)'),2)+' minutes'
    ;
    printf,tfile,'Time for Base & Sat Fix was '+strtrim(string(float(T_first-T_start)/60.0,format='(f12.3)'),2)+' minutes'
    printf,tfile,'Time for Filter was '+strtrim(string(float(T_sec-T_first)/60.0,format='(f12.3)'),2)+' minutes'
    printf,tfile,'Time for Base Re-Fix was '+strtrim(string(float(T_third-T_sec)/60.0,format='(f12.3)'),2)+' minutes'
    printf,tfile,'Time for anc2at was '+strtrim(string(float(T_end-T_third)/60.0,format='(f12.3)'),2)+' minutes'
    printf,tfile,'Total Elapsed Time was '+strtrim(string(float(T_end-T_start)/60.0,format='(f12.3)'),2)+' minutes'
    flush,tfile
    
    ;check open files
    fids=envi_get_file_ids()
    if(fids[0] ge 0) then begin
      for i=0,n_elements(fids)-1 do begin
        envi_file_query,fids[i],fname=tname
        print,'Open File '+strtrim(string(i+1),2)+'='+strtrim(tname,2)
        envi_file_mng,id=fids[i],/remove
      endfor
    endif
    
    help,/memory
    
  endfor
  
  T_total=systime(1)
  printf,tfile,'TOTAL Elapsed Time for all Cases was '+strtrim(string(float(T_total-T_outer)/60.0,format='(f12.3)'),2)+' minutes'
  flush,tfile
  
  free_lun,tfile,/force
  logfile_set=0b
  
  out:
  ;
  free_lun,tfile,/force
  logfile_set=0b
  ;
  print,'Nu_Script_Batch_pulse_filt finished'
  if (err_flag) then begin
    print,'An ERROR has occurred! Check the Log File'
    print,'Log File Name='+strtrim(logfile,2)
  endif
  heap_gc,/verbose
  if (batch_mode and (err_flag eq 0) and (bale_out eq 0)) then begin
    envi_batch_exit,/exit_idl,/no_confirm
  endif
  
end
