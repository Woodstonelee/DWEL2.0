; a script to do a batch job of writing full AT Projected files
pro Nu_Script_Batch_test_resamp
  compile_opt idl2
  
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
    goto, out
  endif
  
  print,' '
  print,'Starting Nu_Script_Batch_Pt_Cloud'
  print,' '
  
  envi, /restore_base_save_files
  help,/memory
  
  ;clean up any fids which are no longer where they were!
  ;ENVI issue that is annoying and leads to confusion
  clean_envi_file_fids
  
  ;
  bale_out=0b
  err_flag=0b
  batch_mode=0b
  
  incfile=[$
    'Y:\DWEL\Data\Processed_data\TumEE05_Ptcld\TumEE05_waveform_2014-08-07-13-37_1548_cube_pulse_filt_update.img',$
    'Y:\DWEL\Data\Processed_data\TumEE05_Ptcld\TumEE05_waveform_2014-08-07-13-37_1064_cube_pulse_filt_update.img' $
    ]
  ancfile=[$
    'Y:\DWEL\Data\Processed_data\TumEE05_Ptcld\TumEE05_waveform_2014-08-07-13-37_1548_cube_pulse_filt_update_ancillary.img',$
    'Y:\DWEL\Data\Processed_data\TumEE05_Ptcld\TumEE05_waveform_2014-08-07-13-37_1064_cube_pulse_filt_update_ancillary.img' $
    ]
  ;outpath is the path for the output files - names set up automatically
  outpath=[$
    'Y:\DWEL\Data\Processed_data\TumEE05_resamp_test\',$
    'Y:\DWEL\Data\Processed_data\TumEE05_resamp_test\' $
    ;'Y:\DWEL\Data\Processed_data\Pye2_Ptcld_test\',$
    ;'Y:\DWEL\Data\Processed_data\Pye2_Ptcld_test\' $
    ]
    
  zen_tweak=[$
    -400L,-400L $
    ;3154,3154 $ ;This is pye2
    ]
    
  ;point cloud settings
  runcode=0
  cal_dat=0b
  dwel_az_n=0.0
  save_zero_hits=1
  add_dwel=0
  zhigh=50.0
  zlow=-5.0
  xmin=-50.0
  xmax=50.0
  ymin=-50.0
  ymax=50.0
  save_br=1b
  
  ncase=n_elements(incfile)
  if ((ncase le 0) or (n_elements(ancfile) ne ncase) or (n_elements(outpath) ne ncase)) then begin
    print,'Number of cases does not match among input files etc!'
    goto,out
  endif
  
  print,'ncase='+strtrim(string(ncase),2)
  
  for k=0,ncase-1 do begin
    incfile[k]=strtrim(incfile[k],2)
    if (~file_test(incfile[k])) then begin
      print,'Input dwel cube file does NOT exist!'
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
  endfor
  
  print,'input checking complete'
  
  ;envi_batch_init if batch_mode is set to 1 (true)
  if (batch_mode) then begin
    envi_batch_init,/NO_STATUS_WINDOW
    envi_batch_status_window,/off
  endif
  
  T_outer=systime(1)
  
  for k=0,ncase-1 do begin
  
    print,''
    print,'Outer Run Case='+strtrim(string(k+1),2)
    print,'Input Cube file='+strtrim(incfile[k],2)
    
    inbase=file_basename(incfile[k])
    nbuf=strpos(inbase,'.img')
    baseout=strtrim(strmid(inbase,0,nbuf),2)
    print,'baseout='+baseout
    ;
    outdir=strtrim(outpath[k],2)
    incube=strtrim(incfile[k],2)
    ;outcube=outdir+baseout+'_ptcl.txt'
    outcube=outdir+baseout+'_test_tzero.img'
    
    print,''
    print,'outcube='+outcube
    
    ;Create the output directory (if it already exists, this code will do nothing)
    o_dir = file_dirname(outcube)
    file_mkdir, o_dir
    
    ;check it is not open in envi
    if(file_test(outcube)) then begin
      fids=envi_get_file_ids()
      if(fids[0] eq -1) then begin
      ;
      endif else begin
        for i=0,n_elements(fids)-1 do begin
          envi_file_query,fids[i],fname=tname
          if (strtrim(strlowcase(outcube),2) eq $
            strtrim(strlowcase(tname),2)) then begin
            envi_file_mng,id=fids[i],/remove
            print,'outcube fid removed'
          endif
        endfor
      endelse
    endif
    
    T_start=systime(1)
    
    ;now the point cloud!
    
    settings={ $
      runcode:runcode,$
      cal_dat:cal_dat,$
      DWEL_az_n:dwel_az_n,$
      save_zero_hits:save_zero_hits,$
      add_dwel:add_dwel,$
      save_br:save_br,$
      zhigh:zhigh,$
      zlow:zlow,$
      xmin:xmin,$
      xmax:xmax,$
      ymin:ymin,$
      ymax:ymax $
      }
      
    ;help,settings
    ;print,settings
      
    err=0
    get_info_stats=1b
    
    ;DWEL_get_point_cloud,incfile[k],ancfile[k],outcube,settings,err
    DWEL_Filtered_FixBase_Cmd_oz, incfile[k],ancfile[k],outcube, get_info_stats, zen_tweak[k], err
    
    if (err eq 0) then begin
      print,'point cloud finished without error set'
    endif else begin
      print,'point cloud finished with error set'
      goto,out
    endelse
    
    ;outname1=strtrim(strmid(outcube,0,strpos(outcube, '.', /reverse_search))+'_points.txt',2)
    ;outname2=strtrim(strmid(outcube,0,strpos(outcube, '.', /reverse_search))+'_metadata.txt',2)
    
    ;if(~file_test(outname1)) then begin
    ;  print,'Output points file from DWEL_get_point_cloud does NOT exist!'
    ;  print,'expected file is='+strtrim(outname1,2)
    ;  goto,out
    ;endif
    
    ;if(~file_test(outname2)) then begin
    ;  print,'Output metadata file from DWEL_get_point_cloud does NOT exist!'
    ;  print,'expected file is='+strtrim(outname2,2)
    ;  goto,out
    ;endif
    
    T_end=systime(1)
    print,' '
    print,'Total Elapsed Time was '+strtrim(string(float(T_end-T_start)/60.0,format='(f12.3)'),2)+' minutes'
    
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
  print,'TOTAL Elapsed Time for all Cases was '+strtrim(string(float(T_total-T_outer)/60.0,format='(f12.3)'),2)+' minutes'
  
  out:
  
  ;
  print,'Nu_Script_Batch_test_resamp finished'
  heap_gc,/verbose
  
  if (batch_mode and (err_flag eq 0) and (bale_out eq 0)) then begin
    envi_batch_exit,/exit_idl,/no_confirm
  endif
  
end
