; a script to do a batch job of writing full AT Projected files, a point cloud
; and a post-point cloud "pfilter" image
pro nu_script_batch_atp_ptcld
compile_opt idl2
;
bale_out=0b
batch_mode=0b
err_flag=0b
logfile_set=0b
tfile=30

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
print,'Starting nu_script_batch_atp_ptcld'
print,' '

envi, /restore_base_save_files
help,/memory

;clean up any fids which are no longer where they were!
;ENVI issue that has been annoying and leads to confusion
clean_envi_file_fids
;
err_flag=0b
err=0

;====================================================================
;settings section for nu_script_batch_atp_ptcld
;set up the name for the log file
logfile='Y:\DWEL\Data\Processed_data\wire10092014\pulse_atp_ptcld_wire10092014_test_logfile.log'
;
;incfile is the file to be converted to AT project and a point cloud obtained
incfile=[$
'Y:\DWEL\Data\Processed_data\wire10092014\wire10_waveform_2014-09-10-13-44_1548_cube_pulse_filt_update.img',$
'Y:\DWEL\Data\Processed_data\wire10092014\wire10_waveform_2014-09-10-13-44_1064_cube_pulse_filt_update.img',$
'Y:\DWEL\Data\Processed_data\BE004_03092014_Ptcl_Test\BE004_waveform_2014-09-03-13-37_1548_cube_pulse_filt_update.img',$
'Y:\DWEL\Data\Processed_data\BE004_03092014_Ptcl_Test\BE004_waveform_2014-09-03-13-37_1064_cube_pulse_filt_update.img',$
'Y:\DWEL\Data\Processed_data\TumEE05_07082014\TumEE05_waveform_2014-08-07-13-37_1548_cube_pulse_filt_update.img',$
'Y:\DWEL\Data\Processed_data\TumEE05_07082014\TumEE05_waveform_2014-08-07-13-37_1064_cube_pulse_filt_update.img',$
'Y:\DWEL\Data\Processed_data\Pye2_Ptcld\pye2a_waveform_2014-06-21-10-55_1548_cube_pulse_filt_update.img',$
'Y:\DWEL\Data\Processed_data\Pye2_Ptcld\pye2a_waveform_2014-06-21-10-55_1064_cube_pulse_filt_update.img' $
]
ancfile=[$
'Y:\DWEL\Data\Processed_data\wire10092014\wire10_waveform_2014-09-10-13-44_1548_cube_pulse_filt_update_ancillary.img',$
'Y:\DWEL\Data\Processed_data\wire10092014\wire10_waveform_2014-09-10-13-44_1064_cube_pulse_filt_update_ancillary.img',$
'Y:\DWEL\Data\Processed_data\BE004_03092014_Ptcl_Test\BE004_waveform_2014-09-03-13-37_1548_cube_pulse_filt_update_ancillary.img',$
'Y:\DWEL\Data\Processed_data\BE004_03092014_Ptcl_Test\BE004_waveform_2014-09-03-13-37_1064_cube_pulse_filt_update_ancillary.img',$
'Y:\DWEL\Data\Processed_data\TumEE05_07082014\TumEE05_waveform_2014-08-07-13-37_1548_cube_pulse_filt_update_ancillary.img',$
'Y:\DWEL\Data\Processed_data\TumEE05_07082014\TumEE05_waveform_2014-08-07-13-37_1064_cube_pulse_filt_update_ancillary.img',$
'Y:\DWEL\Data\Processed_data\Pye2_Ptcld\pye2a_waveform_2014-06-21-10-55_1548_cube_pulse_filt_update_ancillary.img',$
'Y:\DWEL\Data\Processed_data\Pye2_Ptcld\pye2a_waveform_2014-06-21-10-55_1064_cube_pulse_filt_update_ancillary.img' $
]
;outpath is the path for the output files - names set up automatically
outpath=[$
'Y:\DWEL\Data\Processed_data\wire10092014\',$
'Y:\DWEL\Data\Processed_data\wire10092014\',$
'Y:\DWEL\Data\Processed_data\BE004_03092014_Ptcl_Test\',$
'Y:\DWEL\Data\Processed_data\BE004_03092014_Ptcl_Test\',$
'Y:\DWEL\Data\Processed_data\TumEE05_07082014\',$
'Y:\DWEL\Data\Processed_data\TumEE05_07082014\',$
'Y:\DWEL\Data\Processed_data\Pye2_Ptcld\',$
'Y:\DWEL\Data\Processed_data\Pye2_Ptcld\' $
]
;zen_tweak corrects for variations in "top" point in scans
zen_tweak=[$
-1300L,-1300L $ ;saturated CLab scan
,-700L,-700L $ ;Vinyard BE004
,-400L,-400L $ ;TumEE04
 ,3154,3154 $  ;Pye2
]
;dwel_az_n is the dwel azimuth corresponding to North
dwel_az_n=[$
0.0,0.0 $
,0.0,0.0 $
,0.0,0.0 $
,0.0,0.0 $
]
;settings for flow control
bale_out=0b
batch_mode=1b
;settings for processes
;at project full
Max_Zenith_Angle=117.0
output_resolution=4.0
;point cloud settings
runcode=0     ;set to a value to distinguish runs
cal_dat=0b    ;cal_dat=1b if you have calibration for the DWEL
save_zero_hits=1 ;if 1 then no hits is recorded as a valid gap
add_dwel=0    ;if set 1 (0,0,0) and (0,0,dwel_height) are set for reference
;set a box of limits for impossible or unnecessary points
;useful to remove impossible points
zhigh=50.0
zlow=-5.0
xmin=-50.0
xmax=50.0
ymin=-50.0
ymax=50.0
;if save_br is set (=1b) the B and R images are saved
;Note these are real valued so very large - used only for debugging
save_br=0b
;save_pfilt=1 saves the Pfilter image as well - normally 1b
save_pfilt=1b

;====================================================================
;
;set up the log file
logfile_set=0b
if (strlen(logfile) le 0) then begin
  print,'the logfile name is empty!!'
  istat=1
  err_flag=1b
  goto,out
endif
;make sure path to logfile exists (can be used to set it up)
o_dir = file_dirname(logfile)
file_mkdir, o_dir

;start basic checking of the input data files
incfile=strtrim(incfile,2)
ncase=n_elements(incfile)
if ((ncase le 0) or (n_elements(ancfile) ne ncase) or (n_elements(outpath) ne ncase) $
    or (n_elements(zen_tweak) ne ncase) or (n_elements(dwel_az_n) ne ncase)) then begin
  print,'Number of cases does not match among input files etc!'
  err_flag=1b
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
  if (DWEL_headers.headers_present le 0s or ~DWEL_headers.run_present) then begin
    print,strtrim('Input file is NOT a DWEL file',2)
    print,'Case='+strtrim(string(k),2)
    print,'File='+strtrim(incfile[k],2)
    envi_file_mng,id=infile_fid,/remove
    err_flag=1b
    goto,out
  endif
;
  if (~DWEL_headers.filtfix_present) then begin
    print,strtrim('Input file is NOT filtered and refixed!',2)
    print,'Case='+strtrim(string(k),2)
    print,'File='+strtrim(incfile[k],2)
    envi_file_mng,id=infile_fid,/remove
    err_flag=1b
    goto,out
  endif
  envi_file_mng,id=infile_fid,/remove
endfor

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
printf,tfile,strtrim('Running at project, point cloud and saving pfilter image for DWEL processing',2)
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
;now list the zen_tweaks
buf=strtrim(string(zen_tweak),2)
buf='['+strtrim(strjoin(buf,','),2)+']'
printf,tfile,'Input Zenith Tweaks='+buf

;print,buf
buf=''

;all has been set up and checked as far as possible
printf,tfile,'All Set and Ready to run for '+strtrim(string(ncase),2)+' runs'
flush,tfile

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
print,'Input Ancillary file='+strtrim(ancfile[k],2)

printf,tfile,'Outer Run Case='+strtrim(string(k+1),2)
printf,tfile,'Input Cube file='+strtrim(incfile[k],2)
printf,tfile,'Input Ancillary file='+strtrim(ancfile[k],2)
flush,tfile

inbase=file_basename(incfile[k])
nbuf=strpos(inbase,'.img')
baseout=strtrim(strmid(inbase,0,nbuf),2)
print,'baseout='+baseout
;
outdir=strtrim(outpath[k],2)
incube=strtrim(incfile[k],2)
outcube=outdir+baseout+'_at_proj_full.img'
outcloud=outdir+baseout+'_atp_ptcl.txt'


print,''
print,'outcube='+outcube
print,'outcloud='+outcloud

; Create the output directory (if it already exists, this code will do nothing)
  o_dir = file_dirname(outcube)
  file_mkdir, o_dir
  o_dir = file_dirname(outcloud)
  file_mkdir, o_dir

;check they are not open in envi
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
;check they are not open in envi
  if(file_test(outcloud)) then begin
    fids=envi_get_file_ids()
    if(fids[0] eq -1) then begin
;
    endif else begin
      for i=0,n_elements(fids)-1 do begin
        envi_file_query,fids[i],fname=tname
        if (strtrim(strlowcase(outcloud),2) eq $
            strtrim(strlowcase(tname),2)) then begin
            envi_file_mng,id=fids[i],/remove
            print,'outcloud fid removed'
        endif
      endfor
    endelse
  endif
  
T_start=systime(1)

printf,tfile,'Running AT Project full procedure'
flush,tfile

err=0
;run the at project
dwel_cube2at_oz,incube,ancfile[k],outcube,Max_Zenith_Angle,output_resolution,zen_tweak[k],err

if (err gt 0) then begin
  print,'DWEL_cube2at_oz returned with error'
  print,'Local error number='+strtrim(string(err),2)
  printf,tfile,'DWEL_cube2at_oz returned with error'
  printf,tfile,'Case Number='+strtrim(string(k+1),2)
  printf,tfile,'Local error number='+strtrim(string(err),2)
  flush,tfile
  err_flag=1b
  goto,out
endif

;test that it worked in practice!
if(~file_test(outcube)) then begin
  print,'Output file from DWEL_cube2at_oz does NOT exist!'
  print,'expected file is='+strtrim(outcube,2)
  printf,tfile,'Output file from DWEL_cube2at_oz does NOT exist!'
  printf,tfile,'Case Number='+strtrim(string(k+1),2)
  printf,tfile,'expected file is='+strtrim(outcube,2)
  flush,tfile
  err_flag=1b
  goto,out
endif

AncillaryFile = strtrim(strmid(outcube,0,strpos(outcube, '.', /reverse_search))+'_extrainfo.img',2)
if(~file_test(AncillaryFile)) then begin
  print,'Expected output Ancillary file does NOT exist!'
  print,'expected file is='+strtrim(AncillaryFile,2)
  printf,tfile,'Expected output Ancillary file does NOT exist!'
  printf,tfile,'Case Number='+strtrim(string(k+1),2)
  printf,tfile,'expected file is='+strtrim(AncillaryFile,2)
  flush,tfile
  err_flag=1b
  goto,out
endif

printf,tfile,'Output cube2at Image file='+strtrim(outcube,2)
printf,tfile,'Output cube2at Ancillary Image='+strtrim(AncillaryFile,2)
flush,tfile

;now the point cloud!
;settings are controlled at the beginning of this file
settings={ $
runcode:runcode+k,$
cal_dat:cal_dat,$
DWEL_az_n:dwel_az_n[k],$
save_zero_hits:save_zero_hits,$
add_dwel:add_dwel,$
save_br:save_br,$
save_pfilt:save_pfilt,$
zhigh:zhigh,$
zlow:zlow,$
xmin:xmin,$
xmax:xmax,$
ymin:ymin,$
ymax:ymax $
}

printf,tfile,'Running Point Cloud procedure'
flush,tfile

err=0
;run point cloud routine
DWEL_get_point_cloud,outcube,AncillaryFile,outcloud,settings,err

if (err eq 0) then begin
  print,'point cloud finished without error set'
endif else begin
  print,'point cloud finished with error set'
  printf,tfile,'point cloud finished with error set'
  printf,tfile,'Case Number='+strtrim(string(k+1),2)
  flush,tfile
  err_flag=1b
  goto,out
endelse

outname1=strtrim(strmid(outcloud,0,strpos(outcloud,'.',/reverse_search))+'_points.txt',2)
outname2=strtrim(strmid(outcloud,0,strpos(outcloud,'.',/reverse_search))+'_metadata.txt',2)

if(~file_test(outname1)) then begin
  print,'Output points file from DWEL_get_point_cloud does NOT exist!'
  print,'expected file is='+strtrim(outname1,2)
  printf,tfile,'Output points file from DWEL_get_point_cloud does NOT exist!'
  printf,tfile,'Case Number='+strtrim(string(k+1),2)
  printf,tfile,'expected file is='+strtrim(outname1,2)
  flush,tfile
  err_flag=1b
  goto,out
endif

if(~file_test(outname2)) then begin
  print,'Output metadata file from DWEL_get_point_cloud does NOT exist!'
  print,'expected file is='+strtrim(outname2,2)
  printf,tfile,'Output metadata file from DWEL_get_point_cloud does NOT exist!'
  printf,tfile,'Case Number='+strtrim(string(k+1),2)
  printf,tfile,'expected file is='+strtrim(outname2,2)
  flush,tfile
  err_flag=1b
  goto,out
endif

if (save_pfilt) then begin
  outname3=strtrim(strmid(outcloud,0,strpos(outcloud,'.',/reverse_search))+'_pfilter.img',2)
  if(~file_test(outname3)) then begin
    print,'Output Pfilter file from DWEL_get_point_cloud does NOT exist!'
    print,'expected file is='+strtrim(outname3,2)
    printf,tfile,'Output Pfilter file from DWEL_get_point_cloud does NOT exist!'
    printf,tfile,'Case Number='+strtrim(string(k+1),2)
    printf,tfile,'expected file is='+strtrim(outname3,2)
    flush,tfile
    err_flag=1b
    goto,out
  endif
endif

outname4=strtrim(strmid(outcloud,0,strpos(outcloud,'.',/reverse_search))+'_pcinfo.img',2)
outname5=strtrim(strmid(outcloud,0,strpos(outcloud,'.',/reverse_search))+'_pcinfo_header_list.txt',2)

if(~file_test(outname4)) then begin
  print,'Output pcinfo file from DWEL_get_point_cloud does NOT exist!'
  print,'expected file is='+strtrim(outname4,2)
  printf,tfile,'Output Pcinfo file from DWEL_get_point_cloud does NOT exist!'
  printf,tfile,'Case Number='+strtrim(string(k+1),2)
  printf,tfile,'expected file is='+strtrim(outname4,2)
  flush,tfile
  err_flag=1b
  goto,out
endif else begin
  err=0
;now get the headers list from pcinfo as it is now the complete set of dwel processing info
  istat=list_dwel_headers(outname4,outname5,err)
  if (istat or (err gt 0) or ~file_test(outname5)) then begin
    print,'return from list_dwel_headers of pc_info file with error'
    print,'expected file is='+strtrim(outname5,2)
    print,'error number='+strtrim(string(err),2)
    err_flag=1b
    goto,out
  endif
endelse

printf,tfile,'Output Point cloud File='+strtrim(outname1,2)
printf,tfile,'Output Metadata File='+strtrim(outname2,2)
if (save_pfilt) then printf,tfile,'Output Pfilter File='+strtrim(outname3,2)
printf,tfile,'Output Pcinfo File='+strtrim(outname4,2)
printf,tfile,'Output Headers List File='+strtrim(outname5,2)
;
flush,tfile

T_end=systime(1)
print,' '
print,'Total Elapsed Time for Case was '+strtrim(string(float(T_end-T_start)/60.0,format='(f12.3)'),2)+' minutes'
printf,tfile,'Total Elapsed Time for Case was '+strtrim(string(float(T_end-T_start)/60.0,format='(f12.3)'),2)+' minutes'
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
print,'TOTAL Elapsed Time for all Cases was '+strtrim(string(float(T_total-T_outer)/60.0,format='(f12.3)'),2)+' minutes'
printf,tfile,'TOTAL Elapsed Time for all Cases was '+strtrim(string(float(T_total-T_outer)/60.0,format='(f12.3)'),2)+' minutes'
flush,tfile
free_lun,tfile,/force
logfile_set=0b

out:
;
print,'nu_script_batch_atp_ptcld finished'

free_lun,tfile,/force
logfile_set=0b
heap_gc,/verbose
if (batch_mode and (err_flag eq 0) and (bale_out eq 0)) then begin
  envi_batch_exit,/exit_idl,/no_confirm
endif

end
