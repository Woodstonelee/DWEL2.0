pro Test_bits
compile_opt idl2

inlun=30
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
print,'Starting test page'
print,' '

envi, /restore_base_save_files
help,/memory

;clean up any fids which are no longer where they were!
;ENVI issue that is annoying and leads to confusion
clean_envi_file_fids

;scratch area
goto,config

config_file='Y:\DWEL\Data\Raw_data\Labtest_26.6.14_No_Lasers\config_2010-01-15-22-57-31.cfg'

  if (~file_test(config_file)) then begin
    print,'input config file does NOT exist!'
    print,'File='+strtrim(config_file,2)
    goto,out
  endif

print,'configuration file='+strtrim(config_file,2)

openr,inlun,config_file,/get_lun,error=err
if (err ne 0) then begin
  print,'error opening file'
  print,'File='+strtrim(config_file,2)
  goto,out
endif

print,'file open, inlun='+strtrim(string(inlun),2)
buf=''

nrec=0
reading:
if(eof(inlun)) then goto,done
readf,inlun,buf
buf=strtrim(buf,2)
;print,buf
if(buf eq '') then goto,reading
npl=strpos(buf,' = ')
if (npl gt -1) then begin
  left=strmid(buf,0,npl)
  right=strmid(buf,npl+3)
  buf=strtrim(left,2)+'='+strtrim(right,2)
endif
if (nrec le 0) then con_entries=[buf] else con_entries=[con_entries,buf]
nrec=nrec+1
goto,reading
done:

free_lun,inlun,/force

print,'number of records='+strtrim(string(nrec),2)
print,''

for k=0,nrec-1 do begin
con_entries[k]=strtrim(con_entries[k],2)
print,con_entries[k]
endfor

goto,out

next:

Casing_Range=[170.0,180.0]

buf=''
print,string(casing_range)
help,string(casing_range)
print, 'number of elements=',n_elements(string(casing_range))
buf=strtrim(string(casing_range),2)
;buf=strsplit(buf,/extract)
help,buf
print,buf
buf='['+strtrim(strjoin(buf,','),2)+']'
print,'Casing Range(deg) ='+buf

pulses:

;now get Dwel pulses

wavelength=1548

DWEL_pulse_model_dual, wavelength, i_val, t_val, r_val, p_range, p_time, pulse

help,i_val
help,t_val
help, r_val
help,p_range
help,p_time
help,pulse

num_pulse=n_elements(pulse)

outpulse_file='Y:\DWEL\Data\Processed_data\pulse_model\model_pulse_'+strtrim(string(wavelength),2)+'.txt'
  outpath=file_dirname(outpulse_file)
  file_mkdir,outpath

print,'configuration file='+strtrim(outpulse_file,2)

;first see if a file like this exists and is open in envi
  if (file_test(outpulse_file)) then begin
    fids=envi_get_file_ids()
    if(fids[0] ne -1) then begin
      for i=0,n_elements(fids)-1 do begin
        envi_file_query,fids[i],fname=tname
        if (strtrim(strlowcase(outpulse_file),2) eq $
          strtrim(strlowcase(tname),2)) then begin
          envi_file_mng,id=fids[i],/remove
        endif
      endfor
    endif
  endif

openw,inlun,outpulse_file,/get_lun,error=err
if (err ne 0) then begin
  print,'error opening file'
  print,'File='+strtrim(outpulse_file,2)
  goto,out
endif

print,'file open, inlun='+strtrim(string(inlun),2)

printf,inlun,strtrim('Getting Model Pulse from DWEL',2)
printf,inlun,strtrim('Run made at: '+systime(),2)
printf,inlun,strtrim('Number of cases='+strtrim(string(2),2),2)
printf,inlun,strtrim('Wavelength=',2)+strtrim(string(wavelength),2)
flush,inlun

buf=''
buf=strtrim(string(i_val),2)
buf=strjoin(buf,',',/single)
printf,inlun,'I_VAL=['+buf+']'
buf=''
buf=strtrim(string(t_val),2)
buf=strjoin(buf,',',/single)
printf,inlun,'T_VAL=['+buf+']'
buf=''
buf=strtrim(string(r_val),2)
buf=strjoin(buf,',',/single)
printf,inlun,'R_VAL=['+buf+']'
printf,inlun,'i,p_range[i],p_time[i],pulse[i]'
flush,inlun
;

buf=''
for i=0,num_pulse-1 do begin
  buf=strtrim(string(i),2)+','+strtrim(string(p_range[i]),2)+','+strtrim(string(p_time[i]),2)+','+strtrim(string(pulse[i]),2)
  printf,inlun,buf
  flush,inlun
endfor
free_lun,inlun,/force
print,'Output file='+strtrim(outpulse_file,2)
  out:
free_lun,inlun,/force

config:

config_file='Y:\DWEL\Data\Raw_data\Labtest_26.6.14_No_Lasers\config_2010-01-15-22-57-31.cfg'
;get configuration info and title
Config_Info=['']
consum=''
  istat=dwel_get_config_info(Config_File,Config_Info,consum)
  if (istat gt 0) then begin
    print,'error in call to dwel_get_config_info'
    print,'Configuration file='+strtrim(config_file,2)
    goto,out
  endif

prefix=''
ind=min([strlen(consum),7])
print,'ind=',ind
prefix=strtrim(strcompress(strmid(consum,0,ind)),2)+'_'
print,'prefix=',prefix
ind2=strpos(prefix,' ')
if (ind2 gt 0) then index=strmid(prefix,0,ind2)+'_'+strmid(prefix,ind2+1)
print,prefix


heap_gc,/verbose
  return
  end
  