; a script to do a batch job of importing DWEL HDF5 files to ENVI cube
; DWEL2Cube_Cmd2, DWEL_H5File, DataCube_File, Wavelength, DWEL_Height, DWEL_ND_Filter, nadirelevshift

pro Nu_Script_Batch_itpulse_filt
compile_opt idl2

bale_out=0b
batch_mode=0b
err_flag=1b
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
print,'Starting Nu_Script_Batch_itpulse_filt'
print,' '

  envi, /restore_base_save_files

help,/memory

;clean up any fids which are no longer where they were!
;ENVI issue that is annoying and leads to confusion
clean_envi_file_fids

;set up the log file
logfile_set=0b
;set up logfile
logfile='Y:\DWEL\Data\Processed_data\TumEE05_07082014_pulse_test\test_pulse_filter_logfile.log'
if (strlen(logfile) le 0) then begin
  print,'the logfile name is empty!!'
  istat=1
  err_flag=1b
  goto,out
endif
;make sure path to logfile exists (can be used to set it up)
o_dir = file_dirname(logfile)
file_mkdir, o_dir

;====================================================================
;incfile is the file to be filtered
incfile=[$
'Y:\DWEL\Data\Processed_data\TumEE05_07082014\TumEE05_waveform_2014-08-07-13-37_1064_cube_bsfix_image.img',$
'Y:\DWEL\Data\Processed_data\TumEE05_07082014\TumEE05_waveform_2014-08-07-13-37_1548_cube_bsfix_image.img' $
;'Y:\DWEL\Data\Processed_data\Pye2_Ptcld\pye2a_waveform_2014-06-21-10-55_1064_cube_pulse_filt_update_at_proj_full.img',$
;'Y:\DWEL\Data\Processed_data\Pye2_Ptcld\pye2a_waveform_2014-06-21-10-55_1548_cube_pulse_filt_update_at_proj_full.img' $
]
;outpath is the path for the output files - names set up automatically
outpath=[$
'Y:\DWEL\Data\Processed_data\TumEE05_07082014_pulse\',$
'Y:\DWEL\Data\Processed_data\TumEE05_07082014_pulse\' $
;'Y:\DWEL\Data\Processed_data\Pye2_Ptcld\',$
;'Y:\DWEL\Data\Processed_data\Pye2_Ptcld\' $
]
;inwl is the wavelength of the laser
inwl=[$
1064,1548$;,1064,1548$
]

;settings for flow control

bale_out=0b
batch_mode=0b
err_flag=0b

;====================================================================
incfile=strtrim(incfile,2)
ncase=n_elements(incfile)
if ((ncase le 0) or (n_elements(outpath) ne ncase) or (n_elements(inwl) ne ncase)) then begin
  print,'Number of cases bad in one or more of the input arrays of info!'
  err_flag=1b
  goto,out
endif

print,'ncase='+strtrim(string(ncase),2)

for k=0,ncase-1 do begin
  incfile[k]=strtrim(incfile[k],2)
  if (~file_test(incfile[k])) then begin
    print,'an input dwel cube file does NOT exist!'
    print,'Case='+strtrim(string(k),2)
    print,'File='+strtrim(incfile[k],2)
    err_flag=1b
    goto,out
  endif
  if (strpos(incfile[k],strtrim(string(inwl[k]),2)) lt 0) then begin
    print,'the name does not agree with the wavelength!'
    print,'Case='+strtrim(string(k),2)
    print,'File='+strtrim(incfile[k],2)
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
printf,tfile,strtrim('Running Nu_Script_Batch_itpulse_filt plus anc2at for DWEL itfilter processing',2)
printf,tfile,strtrim('Run made at: '+systime(),2)
printf,tfile,strtrim('Number of cases='+strtrim(string(ncase),2),2)
flush,tfile

logfile_set=1b

;write out input files
printf,tfile,'Input  files:'
for j=0L,ncase-1L do begin
  printf,tfile,strtrim(incfile[j],2)
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

inbase=file_basename(incfile[k])
nbuf=strpos(inbase,'.img')
baseout=strtrim(strmid(inbase,0,nbuf),2)
print,'baseout='+baseout

;
outdir=outpath[k]
incube=strtrim(incfile[k],2)
outcube=outdir+baseout+'_itpulse_filt.img'

print,''
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

wavelength=inwl[k]
if (wavelength eq 1064) then swop_wl=1548 else swop_wl=1064

DWEL_pulse_model_dual_oz, swop_wl, i_val, t_val, r_val, p_range, p_time, itpulse, t_fwhm, r_fwhm

print,'wavelength=',swop_wl
print,'i_val=',i_val
print,'t_val=',t_val
print,'r_val=',r_val
print,'t_fwhm=',t_fwhm
print,'r_fwhm=',r_fwhm

help,itpulse

print,''
print,'Number of values in filter='+strtrim(string(n_elements(itpulse)),2)

itpulse=itpulse/total(itpulse)
ierr=0
maxval=max(itpulse,mpos)

print,'Mpos='+strtrim(string(mpos),2)
istat=dwel_general_filter(incube,itpulse,mpos,outcube,ierr)

if(~file_test(outcube) or (ierr ne 0) or (istat ne 0)) then begin
  print,'Output file from dwel_general_filter does NOT exist or error!'
  print,'expected file is='+strtrim(outcube,2)
  printf,tfile,'Output file from dwel_general_filter does NOT exist or error!'
  printf,tfile,'expected file is='+strtrim(outcube,2)
  err_flag=1b
  goto,out
endif

print,'Output iterated filtered Image file written for '+'='+strtrim(outcube,2)
printf,tfile,'Output iterated filtered Image file written for '+'='+strtrim(outcube,2)

T_end=systime(1)
print,' '
print,'Total Elapsed Time was '+strtrim(string(float(T_end-T_start)/60.0,format='(f12.3)'),2)+' minutes'
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

out:

;
free_lun,tfile,/force
logfile_set=0b
;
print,'Nu_Script_Batch_itpulse_filt finished'
if (err_flag) then begin
  print,'An ERROR has occurred! Check the Log File'
  print,'Log File Name='+strtrim(logfile,2)
endif
heap_gc,/verbose
if (batch_mode and (err_flag eq 0) and (bale_out eq 0)) then begin
  envi_batch_exit,/exit_idl,/no_confirm
endif

end
