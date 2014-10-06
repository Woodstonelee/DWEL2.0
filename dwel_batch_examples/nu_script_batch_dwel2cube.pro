; a script to do a batch job of importing DWEL HDF5 files to ENVI cube
pro Nu_Script_Batch_DWEL2Cube
compile_opt idl2

tfile=30
batch_mode=1b
bale_out=0b
err_flag=0b
logfile_set=0b

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
print,'Starting Nu_Script_Batch_DWEL2Cube'
print,' '

envi, /restore_base_save_files

;clean up any fids which are no longer where they were!
;ENVI issue that is annoying and leads to confusion
clean_envi_file_fids

;err_flag=1b indicates an error has occured during processing
err_flag=0b
;logfile set governs error response as well
logfile_set=0b
AncillaryFile=''

;Settings section
;Set bale_out=1 if you wish to just do the first runs (two bands of one hdf file) to be sure!
;set batch=1 if you wish it to run in batch and close IDL at the end

;====================================================================
bale_out=0b
batch_mode=1b
;set up logfile
logfile='Y:\DWEL\Data\Processed_data\BE004_03092014_Ptcl_Test\dwel2cube_BE004_logfile.log'
if (strlen(logfile) le 0) then begin
  print,'the logfile name is empty!!'
  istat=1
  err_flag=1b
  goto,out
endif
;make sure path to logfile exists (can be used to set it up)
o_dir = file_dirname(logfile)
file_mkdir, o_dir
;inhdf is the list of HDF files input
inhdf=[$
'Y:\DWEL\Data\Raw_data\BE004_03092014\waveform_2014-09-03-13-37.hdf5' $
]
;confile is the list of .cfg (configuration) files
confile=[$
'Y:\DWEL\Data\Raw_data\BE004_03092014\config_2009-12-31-22-26-26.cfg' $
]
;outpath is the path for the output files - names set up automatically
outpath=[$
'Y:\DWEL\Data\Processed_data\BE004_03092014_Ptcl_Test\' $
]
;DWEL height is the height to mirror centre from the ground
DWEL_Height=[$
0.86 $
]
;constant settings at this stage for DWEL2Cube
beam_div=2.5
srate=2.0
;max angle here is 180 to look at ALL the data in final anc2at
Max_Zenith_Angle=180.0
output_resolution=2.5
;====================================================================

;now start to get things checked
ncase=n_elements(inhdf)

;nadirrelevshift is a setting you will not need to set! Replaced by tweak
nadirelevshift=lonarr(ncase)

;print,'Number of cases=',ncase
if ((ncase le 0) or (n_elements(outpath) ne ncase) or (n_elements(DWEL_Height) ne ncase) $
    or (n_elements(nadirelevshift) ne ncase)) then begin
  print,'Number of cases does not match the other arrays!'
  err_flag=1b
  goto,out
endif
 
print,'Number of Cases='+strtrim(string(ncase),2)

for k=0,ncase-1 do begin
  inhdf[k]=strtrim(inhdf[k],2)
  if (~file_test(inhdf[k])) then begin
    print,'input hdf file does NOT exist!'
    print,'Case='+strtrim(string(k),2)
    print,'File='+strtrim(inhdf[k],2)
    err_flag=1b
    goto,out
  endif
  if (~file_test(confile[k])) then begin
    print,'input configuration file does NOT exist!'
    print,'Case='+strtrim(string(k),2)
    print,'File='+strtrim(confile[k],2)
    err_flag=1b
    goto,out
  endif
  outpath[k]=strtrim(outpath[k],2)
  if (strlen(outpath[k]) le 0) then begin
    print,'input path is blank!'
    print,'Case='+strtrim(string(k),2)
    print,'File='+strtrim(inhdf[k],2)
    err_flag=1b
    goto,out
  endif
  file_mkdir, outpath[k]
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
printf,tfile,strtrim('Running dwel2cube plus anc2at for DWEL hdf to cube processing',2)
printf,tfile,strtrim('Run made at: '+systime(),2)
printf,tfile,strtrim('Number of cases='+strtrim(string(ncase),2),2)
flush,tfile

;logfile now tested and set with some things written to it
logfile_set=1b

;write out input HDF files
printf,tfile,'Input HDF  files:'
for j=0L,ncase-1L do begin
  printf,tfile,strtrim(inhdf[j],2)
endfor
;write out input Configuration files
printf,tfile,'Input CFG  files:'
for j=0L,ncase-1L do begin
  printf,tfile,strtrim(confile[j],2)
endfor
;write output file paths
printf,tfile,'Output file paths:'
for j=0L,ncase-1L do begin
  printf,tfile,strtrim(outpath[j],2)
endfor
;now list the dwell heights
buf=strtrim(string(DWEL_height),2)
buf='['+strtrim(strjoin(buf,','),2)+']'
printf,tfile,'Input DWEL heights='+buf

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
printf,tfile,'Input HDF file='+strtrim(inhdf[k],2)

print,'Outer Run Case='+strtrim(string(k+1),2)

;get configuration info and title
Config_Info=['']
consum=''
  istat=dwel_get_config_info(Confile[k],Config_Info,consum)
  if (istat gt 0) then begin
    print,'error in call to dwel_get_config_info at case='+strtrim(string(k),2)
    print,'Configuration file='+strtrim(confile[k],2)
    goto,out
  endif

prefix=''
ind=min([strlen(consum),8])
prefix=strtrim(strcompress(strmid(consum,0,ind)),2)+'_'
ind2=strpos(prefix,' ')
if (ind2 gt 0) then index=strmid(prefix,0,ind2)+'_'+strmid(prefix,ind2+1)
print,'prefix='+strtrim(prefix,2)

inbase=strtrim(file_basename(inhdf[k]),2)
nbuf=strpos(inbase,'hdf5')
baseout=prefix+strtrim(strmid(inbase,0,nbuf-1),2)
print,'baseout='+baseout

;run 1 1548nm

print,'starting 1548nm run'

; first run for 1548nm - hdf to cube
; NOTE - Oz DWEL band names are swopped in HDF file
inwl=1064
outwl=1548
outdir=outpath[k]

;set up the output file names
outcube=outdir+baseout+'_'+strtrim(string(outwl),2)+'_cube.img'
DWEL_AT_File=outdir+baseout+'_'+strtrim(string(outwl),2)+'_cube_Anc2AT_project.img'

printf,tfile,'first case - '+strtrim(string(outwl),2)+'nm'
printf,tfile,'outcube='+outcube
printf,tfile,'DWEL_AT_File='+DWEL_AT_File

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

printf,tfile,'Beam Divergence='+strtrim(string(beam_div),2)+' mrad'
printf,tfile,'Sampling Rate='+strtrim(string(srate),2)+' s/ns'

T_start=systime(1)

print,'Calling DWEL2Cube_cmd'

err=0
DWEL2Cube_cmd_oz,inhdf[k],confile[k],outcube,inwl,outwl,dwel_height[k],beam_div,srate,nadirelevshift[k],err

if (err ne 0) then begin
  printf,tfile,'Error in call to DWEL2Cube_cmd_oz from Nu_Script_Batch_DWEL2Cube'
  printf,tfile,'Local error code='+strtrim(string(err),2)
  err_flag=1b
  goto,out
endif

if(~file_test(outcube)) then begin
  printf,tfile,'Output cube file does NOT exist!'
  printf,tfile,'expected file is='+strtrim(outcube,2)
  err_flag=1b
  goto,out
endif

AncillaryFile = strtrim(strmid(outcube,0,strpos(outcube, '.', /reverse_search))+'_ancillary.img',2)
if(~file_test(AncillaryFile)) then begin
  printf,tfile,'Expected output Ancillary file does NOT exist!'
  printf,tfile,'expected file is='+strtrim(AncillaryFile,2)
  err_flag=1b
  goto,out
endif

printf,tfile,'Output cube file written for '+strtrim(string(outwl),2)+'='+strtrim(outcube,2)
printf,tfile,'Output ancillary file written for '+strtrim(string(outwl),2)+'='+strtrim(AncillaryFile,2)


T_sec=systime(1)

;
; next run is to get anc2at projected qlook to check encoders and quality

; Create the output directory (if it already exists, this code will do nothing)
  o_dir = file_dirname(DWEL_AT_File)
  file_mkdir, o_dir

;check it is not open in envi
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
          print,'at_file fid removed'
      endif
    endfor
  endelse
endif

max_zenith_angle_at=max_zenith_angle
;max_zenith_angle_at=117.0

printf,tfile,'Max Zenith Angle='+strtrim(string(max_zenith_angle_at),2)+' deg'
printf,tfile,'Output Resolution='+strtrim(string(output_resolution),2)+' mrad'

print,'Calling dwel_anc2at'

err=0
zen_tweak=0L
dwel_anc2at_oz, AncillaryFile, DWEL_AT_File, max_zenith_angle_at, output_resolution,zen_tweak,err

if (err gt 0) then begin
  print,'dwel_anc2at returned with error'
  print,'Local error number='+strtrim(string(err),2)
  printf,tfile,'dwel_anc2at returned with error'
  printf,tfile,'Local error number='+strtrim(string(err),2)
  err_flag=1b
  goto,out
endif

if(~file_test(DWEL_AT_File)) then begin
  printf,tfile,'Output file from dwel_cube2at does NOT exist!'
  printf,tfile,'expected file is='+strtrim(DWEL_AT_File,2)
  err_flag=1b
  goto,out
endif

printf,tfile,'Output AT file written for '+strtrim(string(outwl),2)+'='+strtrim(DWEL_AT_File,2)

T_end=systime(1)
print,' '
printf,tfile,'Elapsed Time for Cube run was '+strtrim(string(float(T_sec-T_start)/60.0,format='(f12.3)'),2)+' minutes'
printf,tfile,'Elapsed Time for AT was '+strtrim(string(float(T_end-T_sec)/60.0,format='(f12.3)'),2)+' minutes'
printf,tfile,'Total Elapsed Time was '+strtrim(string(float(T_end-T_start)/60.0,format='(f12.3)'),2)+' minutes'

if (bale_out) then begin
  print,'Baling out after one case for testing
  print,' '
  print,'Elapsed Time for Test run was '+strtrim(string(float(T_end-T_start)/60.0,format='(f12.3)'),2)+' minutes'
  printf,tfile,'Baling out after one case for testing
  printf,tfile,'Elapsed Time for Test run was '+strtrim(string(float(T_end-T_start)/60.0,format='(f12.3)'),2)+' minutes'
  goto,out
endif

;run 2 1064nm
print,''
print,'Starting 1064nm run'
help,/memory

inwl=1548
outwl=1064

outdir=outpath[k]
outcube=outdir+baseout+'_'+strtrim(string(outwl),2)+'_cube.img'
DWEL_AT_File=outdir+baseout+'_'+strtrim(string(outwl),2)+'_cube_Anc2AT_project.img'

print,''
printf,tfile,'second case - '+strtrim(string(outwl),2)+'nm'
printf,tfile,'outcube='+outcube
printf,tfile,'DWEL_AT_File='+DWEL_AT_File

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
        endif
      endfor
    endelse
  endif

T_start=systime(1)

printf,tfile,'Beam Divergence='+strtrim(string(beam_div),2)+' mrad'
printf,tfile,'Sampling Rate='+strtrim(string(srate),2)+' s/ns'

print,'calling DWEL2Cube_cmd'

err=0
DWEL2Cube_cmd_oz,inhdf[k],confile[k],outcube,inwl,outwl,dwel_height[k],beam_div,srate,nadirelevshift[k],err

if (err ne 0) then begin
  printf,tfile,'Error in call to DWEL2Cube_cmd_oz from Nu_Script_Batch_DWEL2Cube'
  printf,tfile,'Local error code='+strtrim(string(err),2)
  err_flag=1b
  goto,out
endif

if(~file_test(outcube)) then begin
  printf,tfile,'Output cube file does NOT exist!'
  printf,tfile,'expected file is='+strtrim(outcube,2)
  err_flag=1b
  goto,out
endif

AncillaryFile = strtrim(strmid(outcube,0,strpos(outcube, '.', /reverse_search))+'_ancillary.img',2)
if(~file_test(AncillaryFile)) then begin
  printf,tfile,'Expected output Ancillary file does NOT exist!'
  printf,tfile,'expected file is='+strtrim(AncillaryFile,2)
  err_flag=1b
  goto,out
endif

printf,tfile,'Output cube file written for '+strtrim(string(outwl),2)+'='+strtrim(outcube,2)
printf,tfile,'Output ancillary file written for '+strtrim(string(outwl),2)+'='+strtrim(AncillaryFile,2)

T_third=systime(1)

; Create the output directory (if it already exists, this code will do nothing)
  o_dir = file_dirname(DWEL_AT_File)
  file_mkdir, o_dir

;check it is not open in envi
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
        endif
      endfor
    endelse
  endif

max_zenith_angle_at=max_zenith_angle
;max_zenith_angle_at=117.0

printf,tfile,'Max Zenith Angle='+strtrim(string(max_zenith_angle_at),2)+' deg'
printf,tfile,'Output Resolution='+strtrim(string(output_resolution),2)+' mrad'

print,'Calling dwel_anc2at_oz'

err=0
zen_tweak=0L
dwel_anc2at_oz, AncillaryFile, DWEL_AT_File, max_zenith_angle_at, output_resolution,zen_tweak,err

if (err gt 0) then begin
  print,'dwel_anc2at returned with error'
  print,'Local error number='+strtrim(string(err),2)
  printf,tfile,'dwel_anc2at returned with error'
  printf,tfile,'Local error number='+strtrim(string(err),2)
  err_flag=1b
  goto,out
endif

if(~file_test(DWEL_AT_File)) then begin
  printf,tfile,'Output file from dwel_cube2at does NOT exist!'
  printf,tfile,'expected file is='+strtrim(DWEL_AT_File,2)
  err_flag=1b
  goto,out
endif

printf,tfile,'Output AT file written for '+strtrim(string(outwl),2)+'='+strtrim(DWEL_AT_File,2)


T_end=systime(1)
print,' '
printf,tfile,'Elapsed Time for Cube run was '+strtrim(string(float(T_sec-T_start)/60.0,format='(f12.3)'),2)+' minutes'
printf,tfile,'Elapsed Time for AT was '+strtrim(string(float(T_end-T_sec)/60.0,format='(f12.3)'),2)+' minutes'
printf,tfile,'Total Elapsed Time was '+strtrim(string(float(T_end-T_start)/60.0,format='(f12.3)'),2)+' minutes'

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
printf,tfile,'Nu_Script_Batch_DWEL2Cube finished'
flush,tfile
free_lun,tfile,/force
logfile_set=0b

out:
free_lun,tfile,/force
logfile_set=0b
;
print,'Nu_Script_Batch_DWEL2Cube finished'
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
