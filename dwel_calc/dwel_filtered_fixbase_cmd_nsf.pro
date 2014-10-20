;
;+
; NAME:
;DWEL_FILTERED_FIXBASE_CMD_NSF
;
; PURPOSE:
;Baseline re-fix (and ancillary file update) for pre-filtered DWEL data.
;
; CATEGORY:
;DWEL waveform processing.
;
; CALLING SEQUENCE:
;dwel_filtered_fixbase_cmd_nsf, FilteredFile, Inancfile, OutUpdatedFile, get_info_stats, zen_tweak, err
;
; INPUTS:
;FilteredFile = the file name of the DWEL cube file that had been base and sat
;fixed and then filtered with dwel_general_filter (convolution with a pulse
;model). 
;Inancfile = ancillary file name of the DWEL cube.
;OutUpdatedFile = 
;get_info_stats = a boolean to ask if the program outputs statistics of the
;casing waveforms. 
;zen_tweak = 
;err = 
;
; OUTPUTS:
;
; SIDE EFFECTS:
;None.
;
; RESTRICTIONS:
;None.
;
; PROCEDURE:
;This routine assumes input is filtered basefixed and satfixed file. Inancfile
;will usually be the ancillary file for the basefix and satfix image (which is
;the one filtered in normal operation). It re-assesses the casing information
;and also re-writes satfix ancillary information and headers after this, the
;pulse has a new model of the Iterated Pulse and Tzero should be the position of
;the peak of the outgoing iterated pulse!
;
; MODIFICATION HISTORY:
;David Jupp, Sept 2014 - Created this routine. 
;Zhan Li, Oct 2014 - Added documentation comments.
;-
pro dwel_filtered_fixbase_cmd_nsf, FilteredFile, Inancfile, OutUpdatedFile, get_info_stats, zen_tweak, err
  ;; FilteredFile: the file name of the DWEL cube file that had been base and sat fixed and then filtered
  ;
  compile_opt idl2
  envi, /restore_base_save_files
  envi_batch_init, /no_status_window
  ;
  lun=99
  inlun=105
  ofile=101
  tfile=98
  ctfile=35
  fname=''
  o_name=''
  ctfile=30
  before_casing=0
  check=0b
  if (check) then begin
    out_of_pulse=200
  endif else begin
    out_of_pulse=400
  endelse
  target_dn=512.0
  err=0
  err_flag=0b
  
  ;; distance from casing (edge of casing) to the true Tzero position
  casing2Tzero = 0.065 ; unit: meters
  ;; the FWHM of outgoing pulse, ns
  outgoing_fwhm = 5.1
  ;; the full width of outgoing pulse where intensity is below 0.01 of
  ;;maximum
  pulse_width_range = 40.0
  print,'pulse_width_range=',pulse_width_range
  ;saturation test
  sat_test=1023L
  
  ;clean up any fids which are no longer where they were!
  ;ENVI issue that is annoying and leads to confusion
  clean_envi_file_fids
  
  print,'entering filtered basefix program'
  
  ;set speed of light metres per nsec /2
  c=0.299792458
  c2=c/2.0
  
  ; Open DWEL cube file
  envi_open_file, FilteredFile, r_fid=infile_fid, /no_realize, $
    /no_interactive_query 
  
  if (infile_fid eq -1) then begin
    print,strtrim('Error opening input file',2)
    print,'Input File: '+strtrim(FilteredFile,2)
    err_flag=1b
    envi_file_mng,id=infile_fid,/remove
    err=2
    goto, cleanup
  endif
  
  envi_file_query, infile_fid, ns=ns, nl=nl, nb=nb, wl=wl, $
    xstart=xstart, ystart=ystart, data_type=type, $
    interleave=ftype, fname=fname, dims=dims
    
  x_range=[dims[1],dims[2]]
  y_range=[dims[3],dims[4]]
  
  ;set the type of file
  ft_nam='Unknown'
  case ftype of
    0: ft_nam='BSQ'
    1: ft_nam='BIL'
    2: ft_nam='BIP'
  endcase
  
  ;get number of bytes in input file data
  nbytes=dt2nb(type)
  
  ;get path and DWEL_file name as separate strings
  f_base=file_basename(fname)
  f_path=file_dirname(fname)
  
  ; Open Ancillary file
  if(~file_test(Inancfile)) then begin
    message_text=[ $
      'Ancillary file does not exist!',$
      'Expected Name='+strtrim(Inancfile,2)$
      ]
    print, message_text
    envi_file_mng,id=infile_fid,/remove
    err_flag=1b
    err=3
    goto, cleanup
  endif
  
  envi_open_file, Inancfile, r_fid=ancillaryfile_fid, $
    /no_realize
  if (ancillaryfile_fid eq -1) then begin
    print,strtrim('Error or No opening ancillary file',2)
    print,'Ancillary File='+strtrim(Inancfile,2)
    err_flag=1b
    envi_file_mng,id=infile_fid,/remove
    envi_file_mng,id=ancillaryfile_fid,/remove
    err=4
    goto, cleanup
  endif
  
  envi_file_query, ancillaryfile_fid, nb=nb_anc, nl=nl_anc, ns=ns_anc, data_type=type_anc
  
  if ((nl_anc ne nl) or (ns_anc ne ns) or (nb_anc lt 3)) then begin
    envi_file_mng,id=ancillaryfile_fid,/remove
    print,strtrim('Ancillary Data File does NOT conform with input DWEL Cube !',2)
    print,'Input File: '+strtrim(FilteredFile,2)
    print,'Ancillary File: '+strtrim(Inancfile,2)
    err_flag=1b
    err=5
    goto, cleanup
  endif
  
  print,'ancillary file type=',type_anc
  
  print,'cube file='+strtrim(FilteredFile,2)
  print,'ancillary file='+strtrim(Inancfile,2)
  
  ;now get the DWEL headers that are present
  ;set up a base structure for the DWEL headers
  DWEL_headers={ $
    f_base:f_base $
    }
    
  ;find all of the DWEL headers in the hdr file as defined by FID
  status=DWEL_get_headers(infile_fid,DWEL_headers)
  
  if (not status) then begin
    print,strtrim('Bad FID in DWEL_get_headers! DWEL Header setup cancelled!',2)
    print,'Input File: '+strtrim(FilteredFile,2)
    err_flag=1b
    err=6
    goto, cleanup
  endif
  
  ;help,dwel_headers,/structures
  
  if ((DWEL_headers.headers_present le 0s) or ~DWEL_headers.run_present) then begin
    print,strtrim('Input file is NOT a valid DWEL Cube file!',2)
    print,'Input File: '+strtrim(FilteredFile,2)
    envi_file_mng,id=infile_fid,/remove
    envi_file_mng,id=ancillaryfile_fid,/remove
    err_flag=1b
    err=7
    goto, cleanup
  endif
  
  if (~DWEL_headers.base_present) then begin
    print,strtrim('Input file has NOT been basefixed!',2)
    print,'Input File: '+strtrim(FilteredFile,2)
    envi_file_mng,id=infile_fid,/remove
    envi_file_mng,id=ancillaryfile_fid,/remove
    err_flag=1b
    err=8
    goto, cleanup
  endif
  
  info = DWEL_headers.DWEL_scan_info
  
  base_info=DWEL_headers.dwel_base_fix_info
  
  ;also scale
  nu_scale=1.0
  
  ;set the default sampling rate
  match = -1
  for i=0,n_elements(info)-1 do if (strmatch(info[i],'*Sampling Rate*')) then match=i
  if (match ge 0) then begin
    text=strtrim(info[match],2)
    print,'text=',text
    k=strpos(text,'=')
    l=strpos(text,'smp/ns')
    print,'extract=',strtrim(strmid(text,k+1,l-k-1),2)
    ;  reads,strtrim(strmid(text,k+1,l-k-1),2),var2
    var2=float(strtrim(strmid(text,k+1,l-k-1),2))
    if (var2 gt 0.0) then begin
      srate=var2
      srate_set=1
    endif else begin
      srate=2.0
      srate_set=0
    endelse
  endif else begin
    srate=2.0
    srate_set = 0b
  endelse
  if (match ge 0) then print,'info match for sampling rate= ',strtrim(info[match],2)
  if (~srate_set) then print,'sampling rate not set from DWEL headers (using default)'
  print,'sampling rate=',srate
  
  time_step=1.0/srate
  
  print,'time step=',time_step
  
  ;find the Tzero from early basefix
  match = -1
  for i=0,n_elements(base_info)-1 do if (strmatch(base_info[i],'*Tzero=*')) then match=i
  if (match ge 0) then begin
    text=strtrim(base_info[match],2)
    print,'text=',text
    k=strpos(text,'=')
    Tzero=float(strtrim(strmid(text,k+1),2))
  endif else begin
    Tzero=0.0
  endelse
  if (match ge 0) then print,'info match for Tzero= ',strtrim(base_info[match],2)
  print,'Old Tzero=',Tzero
  
  print,'initial wl[0], d_wl (meter) = '+strtrim(string(wl[0]),2)+','+strtrim(string(wl[1]-wl[0]),2)
  
  wl=wl/c2+Tzero
  wl=findgen(nb)*time_step
  
  print,'changed wl[0], d_wl (ns) = '+strtrim(string(wl[0]),2)+','+strtrim(string(wl[1]-wl[0]),2)
  
  ;set the low zenith range from early basefix
  match = -1
  for i=0,n_elements(base_info)-1 do if (strmatch(base_info[i],'*Low(deg)=*')) then match=i
  if (match ge 0) then begin
    text=strtrim(base_info[match],2)
    print,'text=',text
    k=strpos(text,'=')
    low=float(strtrim(strmid(text,k+1),2))
  endif else begin
    low=170.0
  endelse
  if (match ge 0) then print,'info match for Low= ',strtrim(base_info[match],2)
  print,'Low=',low
  
  ;set the high zenith range from early basefix
  match = -1
  for i=0,n_elements(base_info)-1 do if (strmatch(base_info[i],'*High(deg)=*')) then match=i
  if (match ge 0) then begin
    text=strtrim(base_info[match],2)
    print,'text=',text
    k=strpos(text,'=')
    high=float(strtrim(strmid(text,k+1),2))
  endif else begin
    high=180.0
  endelse
  if (match ge 0) then print,'info match for High= ',strtrim(base_info[match],2)
  print,'High=',high
  
  DWEL_Adaptation = ENVI_GET_HEADER_VALUE(ancillaryfile_fid, 'DWEL_Adaptation', undefined=undef)
  
  ;now get the DWEL wavelength
  match = -1
  info = DWEL_Adaptation
  for i=0,n_elements(info)-1 do begin
    if (strmatch(info[i],'*Wavelength=*', /fold_case)) then match=i
  endfor
  if match ge 0 then begin
    text=strtrim(info[match],2)
    print,'text=',text
    k=strpos(text,'=')
    print,'extract=',strtrim(strmid(text,k+1,4),2)
    wavelength=fix(strtrim(strmid(text,k+1,4),2))
  endif else begin
    if (strpos(DWEL_headers.f_base, '1064') ne -1) then begin
      wavelength = 1064
    endif
    if (strpos(DWEL_headers.f_base, '1548') ne -1) then begin
      wavelength = 1548
    endif
  endelse
  
  if (wavelength eq 1064) then begin
    target_dn=512.0
  endif else target_dn=509.0
  
  ;input some planes of data from the ancillary file
  anc_data=lonarr(ns,nl,9)
  for j=0,8 do begin
    anc_data[*,*,j]=long(ENVI_GET_DATA(fid=ancillaryfile_fid, dims=[-1L,0,ns-1,0,nl-1], pos=j))
  endfor
  
  ;close up the envi files
  envi_file_mng,id=infile_fid,/remove
  envi_file_mng,id=ancillaryfile_fid,/remove
  
  ;compute the mask from the ancillary file
  m=bytarr(ns,nl)
  mask_all=bytarr(ns,nl)
  m = byte(anc_data[*,*,6])
  if (max(m) gt 1b) then begin
    w = where(m lt 255, nbad)
    ;Create the binary mask
    m = replicate(1b,ns,nl)
    if (nbad gt 0) then m[w] = 0b
  endif
  mask_all=m
  m=0b
  
  scan_encoder=fltarr(ns,nl)
  rotary_encoder=fltarr(ns,nl)
  zeniths=fltarr(ns,nl)
  azimuths=fltarr(ns,nl)
  
  scan_encoder=float(anc_data[*,*,2])
  rotary_encoder=float(anc_data[*,*,3])
  
  ;set up a structure and push it onto the heap
  sav={ $
    Nshots:ns,$
    Nscans:nl,$
    ShotZen:scan_encoder,$
    ShotAzim:rotary_encoder $
    }
    
  ;now put the data on the heap with a pointer
  p_stat=ptr_new(sav,/no_copy)
  ;compute the zeniths and azimuths
  ;
  status = dwel_set_theta_phi_nsf(p_stat,zen_tweak)
  ;put the results into the local arrays
  zeniths=(*p_stat).ShotZen
  azimuths=(*p_stat).ShotAzim
  ptr_free, p_stat
  
  ;--------------------------------------------
  ;now get the output file name  
  output:
  ;check if output file is open in envi
  ;if it's open in envi, them close and remove the file from envi. 
  if(file_test(OutUpdatedFile)) then begin
    fids=envi_get_file_ids()
    if(fids[0] eq -1) then begin
    ;
    endif else begin
      for i=0,n_elements(fids)-1 do begin
        envi_file_query,fids[i],fname=tname
        if (strtrim(strlowcase(OutUpdatedFile),2) eq $
          strtrim(strlowcase(tname),2)) then begin
          envi_file_mng,id=fids[i],/remove
        endif
      endfor
    endelse
  endif
  
  ;open the input file
  err=0
  openr,inlun,FilteredFile,/get_lun,error=err
  if (err gt 0) then begin
    err=7
    goto,cleanup
  endif
  
  ;; get mean pulse and baseline from the given casing area designated
  ;;by the zenith angles.
  n = long(0)
  sum = dblarr(nb)
  sum2 = dblarr(nb)
  temp=dblarr(nb)
  line_scale=fltarr(nl)+1.0
  t_zerol=fltarr(nl)
  if (get_info_stats) then save=fltarr(9,nl)
  count=long(0)
  
  ;set up the pointer for read_binary
  bufrs=long64(nbytes)*long64(nb)*long64(ns)
  pointsz=long64(0)
  ;
  for i=0L,nl-1L do begin
    satmask=long(reform(anc_data[*,i,0]))
    index=where((mask_all[*,i] ne 0) and ((zeniths[*,i] ge low) and (zeniths[*,i] le high)) $
      and (satmask ne 1L), count)
    if (count gt 0L) then begin
      ;      data = envi_get_slice(fid=infile_fid, line=i, /bil)
      pointsz=long64(i)*long64(bufrs)
      data=read_binary(inlun,data_start=pointsz,data_dims=[ns,nb],data_type=type)
      d = double(reform(data[index,*]))
      data=0b
      temp=total(d,1,/double)/double(count)
      temp2=total(d[*,out_of_pulse:nb-1],/double)/(double(count)*double(nb-out_of_pulse))
      ;
      nt=n_elements(temp)
      store=reform(temp[before_casing:nt-1])
      tempmax=max(store,nct)
      ; interpolate peak location
      istat = peak_int([float(nct-1), float(nct), float(nct+1)], store[[nct-1, nct, nct+1]], tzero_loc, value, offset)
      
      mpos=before_casing+nct
      line_scale[i]=target_dn/(tempmax-temp2)
      t_zerol[i]=(float(before_casing)+tzero_loc)*time_step
      if (get_info_stats) then begin
        save[0,i]=float(i)
        save[1,i]=float(count)
        ;      save[2,i]=max(temp,mpos)-temp2
        save[2,i]=tempmax-temp2
        save[3,i]=float(mpos)*time_step
        save[4,i]=float(temp2)
        save[5,i]=float(line_scale[i])
        save[6,i]=float(sqrt(total((temp-temp2)^2,/double)))
        ;      save[7,i]=(wl[1]-wl[0])*save[6,i]/(max(temp)-temp2)
        save[7,i]=(wl[1]-wl[0])*save[6,i]/save[2,i]
        save[8,i]=t_zerol[i]
      endif
      if ((i gt 50L) and (i lt nl-9L)) then begin
        n = long(n) + 1L
        sum = sum + temp
        sum2 = sum2 + total(d^2, 1,/double)/double(count)
      endif
    endif else begin
      line_scale[i]=0.0
      t_zerol[i]=0.0
      if (get_info_stats) then begin
        save[0,i]=float(i)
        save[1,i]=0.0
        save[2,i]=0.0
        save[3,i]=0.0
        save[4,i]=0.0
        save[5,i]=1.0
        save[6,i]=0.0
        save[7,i]=0.0
        save[8,i]=0.0
      endif
    endelse
    index=0b
    data=0b
    d=0b
    temp=0b
    temp2=0b
  endfor
  
  ;make some space
  d=0b
  data=0b
  temp=0b
  index=0b
  d=double(n)
  pulse=float(sum/d)
  sig=float(sqrt(abs(sum2-sum^2/d)/d))
  ;help,sig
  ;make space again
  sum=0b
  sum2=0b
  d=0b
  
  ;pulse and sig are length the number of bands
  mean_base=total(pulse[out_of_pulse:nb-1],/double)/double(nb-out_of_pulse)
  mean_base_sig=total(sig[out_of_pulse:nb-1],/double)/double(nb-out_of_pulse)
  if (abs(mean_base) lt 1.0e-3) then cv_base=0.0 else cv_base=100.0*mean_base_sig/mean_base
  
  pulse=pulse-mean_base
  ;casing_power=sqrt(total(pulse^2,/double))
  casing_power=total(pulse,/double)
  if (abs(casing_power) lt 1.0e-6) then casing_fwhm=0.0 else $
    casing_fwhm=(wl[1]-wl[0])*casing_power/max(pulse)
  ;
  ;=======================================================
  ;if get_info_stats set then write out the saved data
    
  o_path=file_dirname(OutUpdatedFile)
  if(get_info_stats) then begin
    ;strtrim(FilteredFile,2)
    n_base=strlen(OutUpdatedFile)
    n_dot=strpos(OutUpdatedFile,'.',/reverse_search)
    ;
    if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
      clog_file=strtrim(OutUpdatedFile,2)+'_filtfix_trace.log'
    endif else begin
      clog_file=strmid(strtrim(OutUpdatedFile,2),0,n_dot)+'_filtfix_trace.log'
    endelse
    print,'log file name=',clog_file
    ;  clog_file=o_path+path_sep()+clog_file
    ;see if the log file exists & remove if it does!
    if(file_test(clog_file)) then begin
      fids=envi_get_file_ids()
      if(fids[0] eq -1) then begin
      
      endif else begin
        for i=0,n_elements(fids)-1 do begin
          envi_file_query,fids[i],fname=tname
          if (strlowcase(strtrim(strlowcase(clog_file),2)) eq $
            strlowcase(strtrim(strlowcase(tname),2))) then begin
            envi_file_mng,id=fids[i],/remove
          endif
        endfor
        
      endelse
    endif
    ;Open Log file
    text_err=0
    openw, ctfile, clog_file,/get_lun,error=text_err
    if (text_err ne 0) then begin
      print,' '
      print,'error opening the casing trace log file'
      print,'Logfile name='+strtrim(clog_file,2)
      print,'DWEL_Filtered_File_FixBase_Cmd terminating'
      print,' '
      err_flag=1b
      err=9
      goto,cleanup
    endif
    ;
    printf,ctfile,strtrim('DWEL calibration casing trace Log File',2)
    printf,ctfile,strtrim('Run made at: '+systime(),2)
    printf,ctfile,strtrim('Input Cube File: '+strtrim(FilteredFile,2),2)
    flush,ctfile
    ;
    pos_sav=where(reform(save[1,*]) gt 0.0,npos_sav)
    if (npos_sav gt 0) then begin
      casing_mean=mean(reform(save[2,pos_sav]),/double)
      casing_stdev=stddev(reform(save[2,pos_sav]),/double)
      casing_cv=100.0*casing_stdev/casing_mean
      pos_mean=mean(reform(save[3,pos_sav]),/double)
      pos_stdev=stddev(reform(save[3,pos_sav]),/double)
      pos_cv=100.0*pos_stdev/casing_mean
      off_mean=mean(reform(save[4,pos_sav]),/double)
      off_stdev=stddev(reform(save[4,pos_sav]),/double)
      off_cv=100.0*pos_stdev/casing_mean
      scale_mean=mean(line_scale[pos_sav],/double)
      scale_cv=stddev(line_scale[pos_sav],/double)
      scale_cv=100.0*scale_cv/scale_mean
    endif else begin
      casing_mean=0.0
      casing_stdev=0.0
      casing_cv=0.0
      pos_mean=0.0
      pos_stdev=0.0
      pos_cv=0.0
      off_mean=0.0
      off_stdev=0.0
      off_cv=0.0
      scale_mean=1.0
      scale_cv=0.0
    endelse
    printf,ctfile,'Stats,Mean,CV(%)'
    outstring='Base_Stats='+strtrim(string(mean_base,format='(f10.3)'),2)+','+ $
      strtrim(string(cv_base,format='(f10.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    outstring='Casing_Stats='+strtrim(string(casing_mean,format='(f10.3)'),2)+','+ $
      strtrim(string(casing_cv,format='(f10.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    outstring='Pos_Stats='+strtrim(string(pos_mean,format='(f10.3)'),2)+','+ $
      strtrim(string(pos_cv,format='(f10.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    outstring='Offset_Stats='+strtrim(string(off_mean,format='(f10.3)'),2)+','+ $
      strtrim(string(off_cv,format='(f10.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    outstring='Scale_Stats='+strtrim(string(scale_mean,format='(f10.3)'),2)+','+ $
      strtrim(string(scale_cv,format='(f10.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    outstring='Mean Casing_Stdev='+strtrim(string(casing_power,format='(f10.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    outstring='Casing_fwhm='+strtrim(string(casing_fwhm,format='(f10.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    printf,ctfile,'Line_Num,Casing_Num,Casing_Max_Val,Casing_Max_Pos,Offset,Line_Scale,Casing_Power,Casing_FWHM'
    flush, ctfile
    ;
    for i=0L,nl-1L do begin
      outstring=strtrim(string(save[0,i],format='(f10.0)'),2)+','+ $
        strtrim(string(save[1,i],format='(f10.0)'),2)+','+ $
        strtrim(string(save[2,i],format='(f10.3)'),2)+','+ $
        strtrim(string(save[3,i],format='(f10.0)'),2)+','+ $
        strtrim(string(save[4,i],format='(f10.3)'),2)+','+ $
        strtrim(string(save[5,i],format='(f10.3)'),2)+','+ $
        strtrim(string(save[6,i],format='(f10.3)'),2)+','+ $
        strtrim(string(save[7,i],format='(f10.3)'),2)+','+ $
        strtrim(string(save[8,i],format='(f10.3)'),2)
      ;
      printf,ctfile,outstring
    endfor
    ;
    flush, ctfile
    free_lun, ctfile,/force
    save=0b
  endif
  ;=======================================================
  
  ; initial time from current data cube before baseline fix
  ;assume wl is now in ns
  time = wl
  ;  p_time = time
  ;  pulse = sum / double(n)
  ;  sig = sqrt((sum2 / double(n) - pulse^2)*double(n)/double(n-1))
  CasingMeanWfMax = max(reform(pulse[before_casing:n_elements(pulse)-1]), nct)
  Tzero_I=before_casing+nct
  print, 'Initial Tzero before baseline fix = ', time[Tzero_I], ' ns'
  tlow=time[Tzero_I] - 1.5*pulse_width_range
  thigh=time[Tzero_I] + 1.5*pulse_width_range
  print,'tlow,thigh=',tlow,thigh
  baseline = dblarr(nb)
  tmpind = where(time lt tlow, count)
  if count gt 0 then begin
    baseline[tmpind] = pulse[tmpind]
  endif
  tmpind = where(time gt thigh, count)
  if count gt 0 then begin
    baseline[tmpind] = pulse[tmpind]
  endif
  tmpind=where((time ge tlow) and (time le thigh),count)
  if (count gt 0) then begin
    baseline[tmpind]=0.0
  endif
  
  pulse = pulse - baseline
  
  ;  CasingMeanWfMax = max(pulse, Tzero_I)
  CasingMeanSig=sig[Tzero_I]
  CasingMeanCV=100.0*CasingMeanSig/CasingMeanWfMax
  
  casing_power=total(reform(pulse[tmpind]),/double)
  tmpmax = max(pulse, Tzero_I)
  if (abs(casing_power) lt 1.0e-6) then casing_fwhm=0.0 else $
    casing_fwhm=(wl[1]-wl[0])*casing_power/max(pulse)
  tmpind=0b
  Tzero=time[Tzero_I]
  print,'Initial Tzero after baseline fix = ',Tzero,' ns'
  
  ;; interpolate peak location
  istat = peak_int(time[[Tzero_I-1, Tzero_I, Tzero_I+1]], pulse[[Tzero_I-1, Tzero_I, Tzero_I+1]], time_int, pulse_int, offset)
  Tzero = float(time_int)
  print, 'Initial Tzero after interpolation = ', Tzero, ' ns'
  
  posscal=where((t_zerol eq 0.0) or (line_scale eq 0.0),nscal)
  if (nscal gt 0) then begin
    t_zerol[posscal]=Tzero
    line_scale[posscal]=scale_mean
  endif
  
  delta= -casing2Tzero/c2 ; 0.065 meter is about the distance between the rotating mirror and the base.
  print,'delta=',delta, ' ns'
  print,''
  ;
  Tzero=Tzero+delta
  print,'Shifted Tzero='+strtrim(string(Tzero),2)+' ns'
  ;now correct the line based estimates as well
  t_zerol=t_zerol+delta
  
  av_tzero=mean(t_zerol[50:nl-11])
  print,'Average Tzero of each scan line = '+strtrim(string(av_tzero),2)+' ns'
  
  t_zerol[0:49]=av_tzero
  t_zerol[nl-10:nl-1]=av_tzero
  
  time=time-Tzero
  
  ;  tend=(Tzero- 1.5*pulse_width_range)>0.0
  ;  tnice=(Tzero- pulse_width_range)>0.0
  ;  shottz=fix(Tzero/time_step)
  ;  shotend=fix(tend/time_step)
  ;  print,'tend=',tend
  ;  print,'shotend,shottz=',shotend,shottz
  
  ;=================================================================================================
  ;write out the information into a spectral library to check later
  write_baseline_file:
  
  ;Determine library file name
  n_base=strlen(OutUpdatedFile)
  n_dot=strpos(OutUpdatedFile,'.',/reverse_search)
  if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
    baseline_library=strtrim(OutUpdatedFile,2)+'_base_Lib.sli'
  endif else begin
    baseline_library=strmid(OutUpdatedFile,0,n_dot)+'_base_Lib.sli'
  endelse
  
  ; check if we need to close and delete library file
  fids=envi_get_file_ids()
  if(fids[0] ne -1) then begin
    for i=0,n_elements(fids)-1 do begin
      envi_file_query,fids[i],fname=name
      if (strtrim(strlowcase(baseline_library),2) eq $
        strtrim(strlowcase(name),2)) then begin
        envi_file_mng,id=fids[i],/remove,/delete
      endif
    endfor
  endif
  
  out_mat=dblarr(nb,3)
  out_mat=double([pulse,baseline,sig])
  out_mat=reform(out_mat,nb,3,/overwrite)
  sense=size(out_mat,/structure)
  
  Lib_Nam=['Pulse','Baseline','sig']
  
  openw, lun, baseline_library,/get_lun
  writeu,lun,out_mat
  free_lun,lun,/force
  
  ;set up the ENVI header for the output library
  descrip='Baseline Fit Information from: '+strtrim(FilteredFile,2)
  envi_setup_head,fname=baseline_library,ns=nb,nl=3,nb=1,$
    xstart=0,ystart=0,file_type=4,interleave=0,$
    data_type=sense.type,spec_names=lib_nam,$
    descrip=descrip, wl=time,$
    zplot_titles=['Time','Value'], $
    /write
    
  ;envi_open_file,baseline_library,r_fid=lib_fid,/no_interactive_query,/no_realize
  ;envi_file_mng,id=lib_fid,/remove
    
  out_mat=0b
  baseline_corr=0b
  ;=================================================================================================
  
  ;set up the DWEL header information for the base fixing
  ;
  dwel_itpulse_model_dual_nsf, wavelength, i_val, t_val, r_val, p_range, p_time, pulse_model
  
  model_fwhm=total(pulse_model)*time_step
  
  DWEL_filtered_fix_info=strarr(21)
  DWEL_filtered_fix_info=[$
    'Program=dwel_filtered_fixbase_cmd_nsf',$
    'Descr=DWEL New Filtered Base Fix Settings',$
    'Processing Date Time='+strtrim(systime(),2),$
    'Pulse='+'NSF_DWEL_Pulse_Model',$
    'Comment=Tzero is the time at which the peak of the output iterated pulse occurs',$
    'Tzero='+strtrim(string(Tzero,format='(f10.3)'),2),$
    'srate='+strtrim(string(srate,format='(f10.2)'),2),$
    'out_of_pulse='+strtrim(string(out_of_pulse,format='(i10)'),2),$
    'Target_dn='+strtrim(string(target_dn,format='(f10.2)'),2),$
    'scale_mean='+strtrim(string(scale_mean,format='(f10.3)'),2),$
    'Noise_RMS='+strtrim(string(mean_base_sig,format='(f10.3)'),2),$
    'Low(deg)='+strtrim(string(low,format='(f10.2)'),2),$
    'High(deg)='+strtrim(string(high,format='(f10.2)'),2),$
    'delta(ns)='+strtrim(string(delta,format='(f10.4)'),2),$
    'casing_max='+strtrim(string(CasingMeanWfMax,format='(f10.3)'),2),$
    'casing_sig='+strtrim(string(CasingMeanSig,format='(f10.3)'),2),$
    'Casing_CV(%)='+strtrim(string(CasingMeanCV,format='(f10.2)'),2),$
    'casing_fwhm(nsec)='+strtrim(string(casing_fwhm,format='(f10.4)'),2),$
    'casing_fwhm(m)='+strtrim(string(casing_fwhm*c2,format='(f10.4)'),2),$
    'model_fwhm(nsec)='+strtrim(string(model_fwhm,format='(f10.4)'),2),$
    'model_fwhm(m)='+strtrim(string(model_fwhm*c2,format='(f10.4)'),2) $
    ]
  DWEL_filtered_fix_info=strtrim(DWEL_filtered_fix_info,2)
  
  ;=====================================================================
  ;now apply the correction to the image
  
  nb_out=nb
  nl_out=nl
  ns_out=ns
  wl_out=fltarr(nb_out)
  wl_range=c2*time
  
  ;Open the output file for BIL tiling
  text_err=0
  ;  openw, ofile, out_name,/get_lun,error=text_err
  ;  if (text_err ne 0) then begin
  ;    print, strtrim('Halting DWEL_baseline_fix', 2)
  ;    print, strtrim(['Error opening output file '+strtrim(out_name,2)], 2)
  ;    goto, cleanup
  ;  endif
  
  ;=====================================================================
  openw, ofile, OutUpdatedFile,/get_lun,error=text_err
  if (text_err ne 0) then begin
    print, strtrim('Halting DWEL_Filtered_FixBase', 2)
    print, strtrim(['Error opening output file '+strtrim(OutUpdatedFile,2)], 2)
    err_flag=1b
    err=10
    goto, cleanup
  endif
  ;=====================================================================
  
  band_pos=indgen(nb_out)
  
  ft_out=1
  ft_str=['BSQ','BIL','BIP']
  
  ;zero=fltarr(ns_out,nb_out)
  temp=fltarr(ns_out,nb_out)
  one_ns=fltarr(ns_out)+1.0
  one_nb=fltarr(nb_out)+1.0
  
  ;  runin_mask=fltarr(nb_out)+1.0
  ;  numpos=indgen(shottz-shotend+1)+shotend
  ;  ord=float(numpos-shotend)/float(shottz-shotend)
  
  ;print,'numpos ends=',numpos[0],numpos[shottz-shotend]
  ;help,ord
  
  ;  print,'numpos'
  ;  print,numpos
  ;  print,''
  ;  print,'ord'
  ;  print,ord
  
  ;  merge=(ord^2)*(3.0-2.0*ord)
  ;  help,merge
  ;  print,'merge'
  ;  print,merge
  ;  runin_mask[0:shotend-1]=0.0
  ;  runin_mask[shotend:shottz]=merge
  ;  runin_mask[shottz+1:nb_out-1]=1.0
  ;  help,runin_mask
  ;  mean_image=fltarr(ns_out,nl_out)
  ;  max_image=fltarr(ns_out,nl_out)
  ;==========================================
  mean_image=fltarr(ns_out,nl_out)
  max_image=fltarr(ns_out,nl_out)
  temp = fltarr(ns_out,nb_out)
  ;==========================================
  ;
  wl_out=fltarr(nb_out)
  posk=where(wl_range lt 0.0,nposk)
  wlk=wl_range[posk[nposk-1]]
  ilambda=-wlk/(time_step*c2)   ; initial fraction of the time between Tzero and
                                ; the closest early waveform bin in a time
                                ; interval
  wl_out=wl_range-wlk
  posk=0b
  posk=where((wl_out ge -6.0) and (wl_out le 95.0),nb_resamp)
  wl_out=reform(wl_out[posk])
  pos_pos=where(wl_out gt 0.0,npos_r)
  if(npos_r le 0) then begin
    print,'Undefined time range, npos_r=',npos_r
    print,'there were NO ranges gt 0.0!'
    err_flag=1b
    err=12
    goto,cleanup
  endif
  fln=1.0/float(npos_r)
  ;
  ;print,'npos_r=',npos_r
  ;help,pos_pos
  
  ;set up the pointer to read the input file again
  bufrs=long(nbytes)*long(nb_out)*long(ns_out)
  pointsz=long64(0)
  
  ;do the processing over the lines for BIL structure
  
  for i=0, nl_out-1 do begin
    ;first get the data tile
    ;    data=envi_get_tile(tile_id,i)
    ;    data=envi_get_slice(fid=infile_fid, line=i, /bil)
    data=read_binary(inlun,data_start=pointsz,data_dims=[ns_out,nb_out],data_type=type)
    pointsz=long64(pointsz)+long64(bufrs)
    pos_z=where(mask_all[*,i] eq 0,count_z)
    temp=float(data)
    
    temp=temp-transpose(baseline+mean_base)##one_ns
    ;    temp=temp*(transpose(runin_mask)##one_ns)
    if (count_z gt 0L) then begin
      temp[pos_z,*]=0.0
    endif
    ;scale the data to standard power level
    temp=scale_mean*temp
    ;
    ;resample to new standard ranges
    wl_loc=fltarr(nb_out)
    wl_loc=wl_range+cmreplicate(c2*(Tzero-t_zerol[i]),nb_out)
    ; ------------------------debug--------------------------
    if (Tzero - t_zerol[i]) gt 0.0 then begin
      if (Tzero - t_zerol[i]) gt ilambda*time_step then begin
        print, 'the bin of zero range is changed from range label at scan line: ' + $
          '', i
      endif 
    endif else begin
      if (t_zerol[i] - Tzero) gt (1-ilambda)*time_step then begin
        print, 'the bin of zero range is changed from range label at scan line: ' + $
          '', i
      endif 
    endelse 
    ; ---------------------end of debug----------------------
    posk=where(wl_loc lt 0.0,nposk)
    wlk=wl_loc[posk[nposk-1]]
    lambda=-wlk/(time_step*c2)
    wl_new=wl_loc-wlk
    posk=0b
    posk=where((wl_new ge -6.0) and (wl_new le 95.0),nb_loc)
    if (nb_loc ne nb_resamp) then begin
      print,'Resampled shot inconsisten, nb_loc='+strtrim(string(nb_loc),2)
      print,'Expected Value='+strtrim(string(nb_resamp),2)
      print,'IDL Line Number='+strtrim(string(i),2)
      print,'Tzero,T_zerol=',Tzero,t_zerol[i]
      err_flag=1b
      err=33
      goto,cleanup
    endif
    ;
    temp=lambda*shift(temp,0,-1)+(1.0-lambda)*temp
    temp=reform(temp[*,posk])
    ;round if integer
    if (type lt 4 or type gt 9) then begin
      temp=fix(round(temp), type=2)
    endif
    ;write out the resulting tile
    writeu,ofile,temp
    ;get some stats images
    mean_image[*,i]=fln*total(temp[*,pos_pos],2)
    max_image[*,i]=float(max(temp[*,pos_pos],DIMENSION=2))
    ;=================================================
    data=0
    temp = 0
    temp=0b
    posk=0b
    wl_loc=0b
  endfor
  
  ;set the output data type
  if (type lt 4 or type gt 9) then begin
    out_type=2
  endif else begin
    out_type=4
  endelse
  
  ;clear up and complete the action
  free_lun, inlun,/force
  free_lun, ofile,/force
  data=0
  temp=0b
  one_ns=0b
  one_nb=0b
  ;  envi_file_mng,id=infile_fid,/remove
  ;===================================================
  ;get names
  out_base=file_basename(OutUpdatedFile)
  ;get output_ancillary file name
  ;Set up the ancillary file
  dot = strpos(OutUpdatedFile,'.',/reverse_search)
  if ((dot lt 0) or ((strlen(OutUpdatedFile)-dot-1) ne 3)) then dot = strlen(OutUpdatedFile)
  ancfile = strmid(OutUpdatedFile, 0, dot)+'_ancillary.img'
  last=strpos(ancfile,path_sep(),/reverse_search)
  anc_base=file_basename(ancfile)
  ;now the maskfile
  dot = strpos(OutUpdatedFile,'.',/reverse_search)
  if ((dot lt 0) or ((strlen(OutUpdatedFile)-dot-1) ne 3)) then dot = strlen(OutUpdatedFile)
  maskfile = strmid(OutUpdatedFile, 0, dot)+'_sat_mask.img'
  mask_base=file_basename(maskfile)
  
  ;===================================================
  ;write out header for the output file
  descrip='DWEL Filter Fix applied to '+strtrim(out_base,2)
  band_names=strarr(nb_resamp)+'Range_Sample_(m)_'
  band_names=band_names+strtrim(string(indgen(nb_resamp)+1),2)+'_filter_fix'
  
  envi_setup_head,fname=OutUpdatedFile,ns=ns_out,nl=nl_out,nb=nb_resamp,$
    xstart=xstart+dims[1],ystart=ystart+dims[3],$
    data_type=out_type, interleave=ft_out, $
    wl=wl_out, inherit=inherit, $
    bnames=band_names,descrip=descrip, $
    zplot_titles=['Range (m)','Intensity'], $
    /write
    
  envi_open_file,OutUpdatedFile,r_fid=out_fid,/no_interactive_query,/no_realize
  
  ;write out the previous header records
  status=DWEL_put_headers(out_fid,DWEL_headers)
  ;
  ;write the new header(s) into the HDR file
  envi_assign_header_value, fid=out_fid, keyword='DWEL_filtered_fix_info', $
    value=DWEL_filtered_fix_info
  envi_write_file_header, out_fid
  
  envi_file_mng,id=out_fid,/remove
  ;===================================================
  print,'ancfile='+strtrim(ancfile,2)
  
  ; check if we need to close and delete ancfile
  fids=envi_get_file_ids()
  if(fids[0] ne -1) then begin
    for i=0,n_elements(fids)-1 do begin
      envi_file_query,fids[i],fname=name
      if (strtrim(strlowcase(ancfile),2) eq $
        strtrim(strlowcase(name),2)) then begin
        envi_file_mng,id=fids[i],/remove,/delete
      endif
    endfor
  endif
  
  ;now write out the new ancillary file
  text_err=0
  openw, ofile, ancfile,/get_lun,error=text_err
  if (text_err ne 0) then begin
    print, strtrim('Halting DWEL_filter_baseline_fix', 2)
    print, strtrim(['Error opening output file '+strtrim(ancname,2)], 2)
    err_flag=1b
    err=13
    goto, cleanup
  endif
  
  pos_pos=where(mask_all ne 0)
  sat_meanmean=mean(mean_image[pos_pos])
  sat_stddevmean=stddev(mean_image[pos_pos])
  sat_mlow=sat_meanmean-4.0*sat_stddevmean
  sat_mhigh=sat_meanmean+4.0*sat_stddevmean
  sat_minmean=min(mean_image[pos_pos]) > sat_mlow
  sat_maxmean=max(mean_image[pos_pos]) < sat_mhigh
  mean_image[pos_pos]=((4095.0*(mean_image[pos_pos]-sat_minmean)/(sat_maxmean-sat_minmean) > 0.0) < 4095.0)
  
  print,'fixed data scale range=',sat_minmean,sat_maxmean
  
  for j=0,4 do begin
    writeu,ofile,long(anc_data[*,*,j])
  endfor
  writeu,ofile,round(max_image)
  ;  writeu,ofile,round(mean_image)
  writeu,ofile,long(mask_all)
  writeu,ofile,round(100.0*zeniths)
  writeu,ofile,round(100.0*azimuths)
  free_lun, ofile,/force
  anc_data=0b
  mean_image=0b
  
  ENVI_SETUP_HEAD, fname=ancfile, $
    ns=ns_out, nl=nl_out, nb=9, $
    interleave=0, data_type=3, $
    /write, $
    bnames=['Sat Mask','Sun Mask','Scan Encoder','Rotary Encoder', $
    'Laser Power','Waveform Mean','Mask','Zenith','Azimuth']
    
  envi_open_file,ancfile,r_fid=anc_fid,/no_interactive_query,/no_realize
  
  ;write out the previous header records
  status=DWEL_put_headers(anc_fid,DWEL_headers)
  ;
  ;write the new header(s) into the HDR file
  envi_assign_header_value, fid=anc_fid, keyword='DWEL_filtered_fix_info', $
    value=DWEL_filtered_fix_info
  envi_write_file_header, anc_fid
  
  envi_file_mng,id=anc_fid,/remove
  
  ;=====================================================================
  ; Do the final cleanup
  
  cleanup:
  free_lun, lun,/force
  free_lun, inlun,/force
  free_lun, ofile,/force
  free_lun, tfile,/force
  free_lun, ctfile,/force
  if (err_flag) then print,'dwel_filtered_fixbase returned with error'
  heap_gc,/verbose
  ;
  return
end