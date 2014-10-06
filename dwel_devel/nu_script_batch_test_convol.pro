; a script to test convol and general filter

pro Nu_Script_Batch_test_convol
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
print,'Starting Nu_Script_Batch_test_convol'
print,' '

  envi, /restore_base_save_files

help,/memory

;clean up any fids which are no longer where they were!
;ENVI issue that is annoying and leads to confusion
clean_envi_file_fids

;
bale_out=0b

incfile=[$
'Y:\DWEL\Data\Processed_data\Test_Final_hdf2cube_speed\PyeSky_waveform_2014-07-02-13-50_1548_cube_bsfix_image.img' $
]

;outpath is the path for the output files - names set up automatically
outpath=[$
'Y:\DWEL\Data\Processed_data\Test_convol\' $
]

inwl=[$
1548 $
]

ncase=n_elements(incfile)
if ((ncase le 0) or (n_elements(outpath) ne ncase) or (n_elements(inwl) ne ncase)) then begin
  print,'Number of cases does not match!'
  goto,out
endif

print,'ncase='+strtrim(string(ncase),2)

for k=0,ncase-1 do begin
  incfile[k]=strtrim(incfile[k],2)
  if (~file_test(incfile[k])) then begin
    print,'an input dwel cube file does NOT exist!'
    print,'Case='+strtrim(string(k),2)
    print,'File='+strtrim(incfile[k],2)
    goto,out
  endif
endfor

print,'input checking complete'
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
outdir=outpath[k]
incube=strtrim(incfile[k],2)
outcube=outdir+baseout+'_convol_filter.img'

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

DWEL_pulse_model_dual_oz,wavelength, i_val, t_val, r_val, p_range, p_time, pulse, t_fwhm, r_fwhm

print,''
print,'Number of values in filter='+strtrim(string(n_elements(pulse)),2)

pulse=pulse/total(pulse)
ierr=0
maxval=max(pulse,mpos)

print,'Mpos='+strtrim(string(mpos),2)

T_mid=systime(1)

err=0
istat=dwel_general_filter(incube,pulse,mpos,outcube,err)

if (err ne 0) then begin
  print,'Error called from dwel_general_filter'
  print,'Error number='+strtrim(string(err),2)
  goto,out
endif

if(~file_test(outcube)) then begin
  print,'Output file from dwel_general_filter does NOT exist!'
  print,'expected file is='+strtrim(outcube,2)
  goto,out
endif

print,'Output Filtered Image file written for '+'='+strtrim(outcube,2)

T_end=systime(1)
print,' '
print,'Total Time for setup was '+strtrim(string(float(T_mid-T_start)/60.0,format='(f12.3)'),2)+' minutes'
print,'Total Time for convol was '+strtrim(string(float(T_end-T_mid)/60.0,format='(f12.3)'),2)+' minutes'
print,'Total Elapsed Time for case was '+strtrim(string(float(T_end-T_start)/60.0,format='(f12.3)'),2)+' minutes'

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
print,'Nu_Script_Batch_test_convol finished'
heap_gc,/verbose
end
