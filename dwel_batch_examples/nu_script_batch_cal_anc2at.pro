; a script to help calibrate zenith point tweak using overlapping areas and anc2at
pro Nu_Script_Batch_cal_anc2at
  compile_opt idl2
  
  bale_out=0b
  err_flag=0b
  err=0
  
  ;Setup protective aunty-Catch to watch out for errors
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
  print,'Starting Nu_Script_Batch_cal_anc2at'
  print,' '
  
  envi, /restore_base_save_files
  
  help,/memory
  
  ;clean up any fids which are no longer where they were!
  ;ENVI issue that is annoying and leads to confusion
  clean_envi_file_fids
  
  err_flag=0b
  err=0
  
  ;
  ;==============================================================================================
  ;nu_script_batch_cal_anc2at settings (no log file in this case)
  bale_out=0b
  
  ;ancfile is the list of ancillary files to use to get anc2at images
  ancfile=[$
    'Y:\DWEL\Data\Processed_data\wire10092014\wire10_waveform_2014-09-10-13-44_1548_cube_ancillary.img',$
    'Y:\DWEL\Data\Processed_data\wire10092014\wire10_waveform_2014-09-10-13-44_1548_cube_ancillary.img',$
    'Y:\DWEL\Data\Processed_data\wire10092014\wire10_waveform_2014-09-10-13-44_1548_cube_ancillary.img' $
    ]
    
  ;outpath is the path for the output files - names set up automatically
  outpath=[$
    'Y:\DWEL\Data\Processed_data\wire10092014\Test_zen_tweak\',$
    'Y:\DWEL\Data\Processed_data\wire10092014\Test_zen_tweak\',$
    'Y:\DWEL\Data\Processed_data\wire10092014\Test_zen_tweak\' $
    ]
    
  zen_tweak=[$
    -1600L,-1300L,-1000L $
    ]
    
  ;============================================================
  ncase=n_elements(ancfile)
  if ((ncase le 0) or (n_elements(outpath) ne ncase) $
    or (n_elements(zen_tweak) ne ncase)) then begin
    print,'Number of cases does not match!'
    err_flag=1b
    goto,out
  endif
  
  print,'ncase='+strtrim(string(ncase),2)
  
  for k=0,ncase-1 do begin
    ancfile[k]=strtrim(ancfile[k],2)
    if (~file_test(ancfile[k]) $
      or ~file_test(ancfile[k],/read)) then begin
      print,'an input dwel cube file does NOT exist!'
      print,'Case='+strtrim(string(k),2)
      print,'Anc File='+strtrim(ancfile[k],2)
      err_flag=1b
      goto,out
    endif
  endfor
  
  print,'input checking complete'
  T_outer=systime(1)
  
  for k=0,ncase-1 do begin
  
    print,''
    print,'Outer Run Case='+strtrim(string(k+1),2)
    print,'Input Cube file='+strtrim(ancfile[k],2)
    
    inbase=file_basename(ancfile[k])
    nbuf=strpos(inbase,'.img')
    baseout=strtrim(strmid(inbase,0,nbuf),2)
    print,'baseout='+baseout
    ;
    outdir=strtrim(outpath[k],2)
    incube=strtrim(ancfile[k],2)
    ;add the run number to the output to distinguish the test images
    outcube=outdir+baseout+'_anc2at_'+strtrim(string(k),2)+'.img'
    
    print,'outcube='+outcube
    
    ; Create the output directory (if it already exists, this code will do nothing)
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
    
    Max_Zenith_Angle=180.0
    output_resolution=3.0
    
    err=0
    ;
    dwel_anc2at_oz, ancfile[k], outcube, max_zenith_angle, output_resolution,zen_tweak[k],err
    ;
    if (err gt 0) then begin
      print,'dwel_cube2at_oz returned with error'
      print,'Local error number='+strtrim(string(err),2)
      err_flag=1b
      goto,out
    endif
    
    if(~file_test(outcube)) then begin
      print,'Output file from dwel_anc2at_oz, does NOT exist!'
      print,'expected file is='+strtrim(outcube,2)
      err_flag=1b
      goto,out
    endif
    
    print,'dwel_anc2at_oz Image file='+strtrim(outcube,2)
    
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
  if (err_flag) then begin
    print,'An ERROR has occurred!'
  endif
  ;
  print,'Nu_Script_Batch_cal_anc2at finished'
  heap_gc,/verbose
end
