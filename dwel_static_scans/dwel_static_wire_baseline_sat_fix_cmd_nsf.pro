;+
; NAME:
;DWEL_STATIC_WIRE_BASELINE_SAT_FIX_CMD_NSF
;
;
; PURPOSE:
;Using wire signal to fix baseline and saturation. Tzero is from wire signals
;shot by shot. 
;
;
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; MODIFICATION HISTORY:
;
;       Mon Nov 3 12:55:43 2014, Zhan Li <zhanli86@geo.bu.edu>
;
;		
;
;-

pro dwel_static_wire_baseline_sat_fix_cmd_nsf, DWELCubeFile, ancillaryfile_name, $
    out_satfix_name, get_info_stats, zen_tweak, err, target_range, settings=settings

  ;; DWELCubeFile: the full file name of the DWEL cube file
    
  compile_opt idl2
;  envi, /restore_base_save_files
;  envi_batch_init, /no_status_window

  print, '************************************************************'
  print, 'The minimum range of stationary scans this processing, '
  print, 'dwel_static_wire_baseline_sat_fix_cmd_nsf can accept = 5.0 m'
  print, 'If your stationary scan is from a closer target, '
  print, 'the processing may give unexpected results'
  print, '************************************************************'

  resolve_routine, 'DWEL_GET_HEADERS', /compile_full_file, /either
  resolve_routine, 'DWEL_SET_THETA_PHI_NSF', /compile_full_file, /either
  resolve_routine, 'DWEL_PULSE_MODEL_DUAL_NSF', /compile_full_file, /either
  resolve_routine, 'DWEL_PUT_HEADERS', /compile_full_file, /either
  resolve_routine, 'CMREPLICATE', /compile_full_file, /either

  print, 'Target range by user: ', target_range

  ;
  lun=99
  inlun=105
  ofile=101
  osatfile=102
  tfile=98
  ctfile=35
  fname=''
  o_name=''
  ctfile=30
  err=0
  
  ;============================================
  ;more internal settings
  ;always get the pulse information and write out files at this time
  get_info_stats=1b
  ;skip these values - currently noisy and useless
  before_casing=100 ; in unit of bins
  ;target dn preset but set later
  target_dn=512.0
  ;used in sun background check to be away from significant return pulses
  far_range=650 ; in unit of bins
  ;threshold on raise of baseline by sun background radiation - note these are scaled and initially basefixed units
  sun_thresh=1.5
  ;the "FWHM" of outgoing pulse (now just a number used and calibrated)
  outgoing_fwhm = 5.1
  ;the full width of outgoing pulse where intensity is below 0.01 of maximum
  pulse_width_range = 40.0 ; in unit of ns
  ;saturation test value in DN
  sat_test=1023L
  ;============================================
  
  finalsettings = { bsfixsettings, $
    ; a priori value of Tzero (scan mirror center, delta corrected), unit, ns,
    ; default is 1064 of NSF DWEL
    ; NSF, 1064, 166.00; 1548, 155.42
    ; Oz, 1064, 225.19; 1548, 218.72
    Tzero:166.00, $
    ; position where you are free of the casing effects - for mean baseline
    ; Oz, 500
    out_of_pulse:400, $
    ; distance between scan mirror center and the wire if present,
    ; unit=meters. Here default for 1064 nm of NSF DWEL.    
    wire2Tzero:0.299}
  ;; if user provides settings
  if n_elements(settings) ne 0 or arg_present(settings) then begin
    finalsettings = update_struct_settings(finalsettings, settings)
  endif 
  out_of_pulse = finalsettings.out_of_pulse
  wire2Tzero = finalsettings.wire2Tzero
  tzero = finalsettings.tzero

  wire2Tzero_user = 0b
  if tag_exist(settings, 'wire2Tzero') then begin
    ;; if user has designate a value for wire2Tzero, then we will stick to
    ;; this value later rather than update it according to instrument and laser
    ;; wavelength. 
    wire2Tzero_user = 1b
  endif 
  tzero_user = 0b
  if tag_exist(settings, 'Tzero') then begin
    ;; if user has designate a value for Tzero, then we will stick to
    ;; this value later rather than update it according to instrument and laser
    ;; wavelength. 
    Tzero_user = 1b
  endif 

  ;; get the size of input file to be processed. It will be used in later
  ;; summary of processing time. 
  procfilesize = file_info(DWELCubeFile)
  procfilesize = procfilesize.size
  ;; get the time now as the start of processing
  starttime = systime(1)

  print,'pulse_width_range='+strtrim(string(pulse_width_range),2)
  
  ;clean up any fids which are no longer where they were!
  ;ENVI issue that is annoying and leads to confusion
  clean_envi_file_fids
  
  print,'entering wire basefix satfix program'
  
  ;set speed of light metres per nsec /2
  c=0.299792458
  c2=c/2.0
  
  ; Open DWEL cube file
  envi_open_file, DWELCubeFile, r_fid=infile_fid,/no_realize,/no_interactive_query
  
  if (infile_fid eq -1) then begin
    print,strtrim('Error opening input file',2)
    print,'Input File: '+strtrim(DWELCubeFile,2)
    err=1
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
  
  ;set number of bytes in input file data
  nbytes=dt2nb(type)
  
  ;get path and DWEL_file name as separate strings
  f_base=file_basename(fname)
  f_path=file_dirname(fname)
  
  ; Open Ancillary file
  if(~file_test(ancillaryfile_name)) then begin
    message_text=[ $
      'Ancillary file is not present',$
      'Expected Name='+strtrim(ancillaryfile_name,2)$
      ]
    print, message_text
    envi_file_mng,id=infile_fid,/remove
    err=2
    goto, cleanup
  endif
  
  envi_open_file, ancillaryfile_name, r_fid=ancillaryfile_fid, $
    /no_realize
  ;check if operation cancelled
  if (ancillaryfile_fid eq -1) then begin
    print,strtrim('Error or No opening ancillary file',2)
    print,'Ancillary File: '+strtrim(ancillaryfile_name,2)
    envi_file_mng,id=infile_fid,/remove
    err=3
    goto, cleanup
  endif
  
  envi_file_query, ancillaryfile_fid, nb=nb_anc, nl=nl_anc, ns=ns_anc, data_type=type_anc
  
  if ((nl_anc ne nl) or (ns_anc ne ns) or (nb_anc lt 3)) then begin
    envi_file_mng,id=ancillaryfile_fid,/remove
    print,strtrim('Ancillary Data File does NOT conform with current DWEL Cube !',2)
    print,'Input File: '+strtrim(DWELCubeFile,2)
    print,'Ancillary File: '+strtrim(ancillaryfile_name,2)
    envi_file_mng,id=infile_fid,/remove
    envi_file_mng,id=ancillaryfile_fid,/remove
    err=4
    goto, cleanup
  endif
  
  print,'ancillary file type=',type_anc
  print,'cube file='+strtrim(DWELCubeFile,2)
  print,'ancillary file='+strtrim(ancillaryfile_name,2)
  
  ;now get the DWEL headers that are present
  ;set up a base structure for the DWEL headers
  DWEL_headers={ $
    f_base:f_base $
    }
    
  ;find all of the DWEL headers in the hdr file as defined by FID
  status=dwel_get_headers(infile_fid,DWEL_headers)
  
  if (~status) then begin
    print,strtrim('Bad FID in DWEL_get_headers! DWEL Header setup cancelled!',2)
    print,'Input File: '+strtrim(DWELCubeFile,2)
    err=5
    goto, cleanup
  endif
  
  ;help,dwel_headers,/structures
  
  if ((DWEL_headers.headers_present le 0s) or (~DWEL_headers.run_present)) then begin
    print,strtrim('Input file is NOT a valid DWEL Cube file!',2)
    print,'Input File: '+strtrim(DWELCubeFile,2)
    envi_file_mng,id=infile_fid,/remove
    envi_file_mng,id=ancillaryfile_fid,/remove
    err=6
    goto, cleanup
  endif
  
  info = DWEL_headers.DWEL_scan_info
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;set baseplate limits
  ;; low=casing_range[0]
  ;; high=casing_range[1]
  low = 0
  high = 0
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ;set the default sampling rate
  match = -1
  for i=0,n_elements(info)-1 do if (strmatch(info[i],'*Sampling Rate*')) then match=i
  if (match ge 0) then begin
    text=strtrim(info[match],2)
    print,'text=',text
    k=strpos(text,'=')
    l=strpos(text,'smp/ns')
    print,'extract=',strtrim(strmid(text,k+1,l-k-1),2)
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
  if (match ge 0) then print,'info match for ',strtrim(info[match],2)
  if (~srate_set) then print,'sampling rate not set'
  print,'sampling rate=',srate
  
  time_step=1.0/srate
  
  print,'time step=',time_step

  ;Read the laser manufacturer from the scan info
  match = -1
  for i=0,n_elements(DWEL_headers.DWEL_scan_info)-1 do begin
    if (strmatch(DWEL_headers.DWEL_scan_info[i],'*lasers*')) then match=i
  endfor
  if (match ge 0) then begin
    sf = strtrim(strcompress(strsplit(DWEL_headers.DWEL_scan_info[match],'=',/extract)),2)
    laser_man = sf[1]
  endif else begin
    laser_man = 'manlight'
  endelse
  
  print,'Laser manufacturer = '+strtrim(laser_man)
  
  ;input some planes of data from the ancillary file
  anc_data=lonarr(ns,nl,9)
  for j=0,8 do begin
    anc_data[*,*,j]=long(envi_get_data(fid=ancillaryfile_fid, dims=[-1L,0,ns-1,0,nl-1], pos=j))
  endfor
  
  DWEL_Adaptation = envi_get_header_value(ancillaryfile_fid, 'DWEL_Adaptation', undefined=undef)
  wavelength = 1064
  if undef then begin
    DWEL_Adaptation = ''
    if (strpos(DWEL_headers.f_base, '1064') ne -1) then begin
      wavelength = 1064
    endif
    if (strpos(DWEL_headers.f_base, '1548') ne -1) then begin
      wavelength = 1548
    endif
    zenithenc = 0L
  endif else begin
    match = -1
    info = DWEL_Adaptation
    for i=0,n_elements(info)-1 do begin
      if (strmatch(info[i],'*Wavelength*', /fold_case)) then match=i
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
    ;; find the scan encoder of zenith point
    for i=0, n_elements(info)-1 do begin
      if (strmatch(info[i], '*Scan encoder of zenith point*', /fold_case)) then $
        match = i
    endfor
    if match ge 0 then begin
      text=strtrim(info[match],2)
      print,'text=',text
      k=strpos(text,'=')
      print,'extract=',strtrim(strmid(text,k+1),2)
      zenithenc=fix(strtrim(strmid(text,k+1),2), type=3)
    endif else begin
      zenithenc = 0L
    endelse 
  endelse

  print, 'wavelength = ' + strtrim(string(wavelength), 2)
  print, 'encoder of zenith point = ' + strtrim(string(zenithenc), 2)

  if (wavelength eq 1064) then begin
    target_dn=512.0
    ;; move leftward from out_of_pulse to find peak of wire signal
    ;;bins_toward_wire = 57 ; in unit of bins
  endif else begin
    target_dn=509.0
    ;; move leftward from out_of_pulse to find peak of wire signal
    ;;bins_toward_wire = 80 ; in unit of bins
  endelse 

  if (~wire2Tzero_user) then begin
    if strcmp(laser_man, 'manlight', /fold_case) then begin
      if wavelength eq 1064 then begin
        wire2Tzero = 0.299 ; unit, meter
      endif else begin
        wire2Tzero = 0.414 ; unit, meter
      endelse 
    endif else begin
      ;; this is Oz DWEL
      if wavelength eq 1064 then begin
        wire2Tzero = 0.156 ; unit, meter
      endif else begin
        wire2Tzero = 0.124 ; unit, meter
      endelse      
    endelse 
  endif 
  print, 'wire2Tzero to be used: ', wire2Tzero, ' m'

  if (~tzero_user) then begin
    if strcmp(laser_man, 'manlight', /fold_case) then begin
      if wavelength eq 1064 then begin
        tzero = 166.00 ; unit, ns
      endif else begin
        tzero = 155.42 ; unit, ns
      endelse 
    endif else begin
      ;; this is Oz DWEL
      if wavelength eq 1064 then begin
        tzero = 225.19 ; unit, ns
      endif else begin
        tzero = 218.72 ; unit, ns
      endelse 
    endelse 
  endif 
  print, 'A priori Tzero to be used: ', tzero, ' ns'
  
  ;close up the envi files
  envi_file_mng,id=infile_fid,/remove
  envi_file_mng,id=ancillaryfile_fid,/remove
  
  ;compute the mask from the ancillary file
  m=bytarr(ns,nl)
  mask_all=bytarr(ns,nl)
  m = byte(anc_data[*,*,6])
  if (max(m) gt 1b) then begin
    print,'max(m) is NOT 1b!!! Fix?'
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
    ShotAzim:rotary_encoder, $
    ZenEnc:zenithenc $
    }
    
  ;now put the data on the heap with a pointer
  p_stat=ptr_new(sav,/no_copy)
  ;compute the zeniths and azimuths
  status = dwel_set_theta_phi_nsf(p_stat,zen_tweak)
  ;put the results into the local arrays
  zeniths=(*p_stat).ShotZen
  azimuths=(*p_stat).ShotAzim
  ptr_free, p_stat
  
  ;------------------------------------------------------------------------------
  output:
  
  ;check if output file is open in envi
  if(file_test(out_satfix_name)) then begin
    fids=envi_get_file_ids()
    if(fids[0] eq -1) then begin
    ;
    endif else begin
      for i=0,n_elements(fids)-1 do begin
        envi_file_query,fids[i],fname=tname
        if (strtrim(strlowcase(out_satfix_name),2) eq $
          strtrim(strlowcase(tname),2)) then begin
          envi_file_mng,id=fids[i],/remove
        endif
      endfor
    endelse
  endif
  
  ;open the input file
  err=0
  ;open input file
  openr,inlun,DWELCubeFile,/get_lun,error=err
  if (err gt 0) then begin
    err=7
    goto,cleanup
  endif
  
  ;; get a mean base from all waveforms
  n = long(0)
  sum = dblarr(nb)
  sum2 = dblarr(nb)
  temp=dblarr(nb)
  line_scale=fltarr(nl)+1.0
  if (get_info_stats) then save=fltarr(8,nl)
  count=long(0)
  ;set up the pointer for read_binary
  bufrs=long(nbytes)*long(nb)*long(ns)
  pointsz=long64(0)
  
  print,'bufrs='+strtrim(string(bufrs),2)
  ;
  ;; find locations of bins for wire signal search
  w_start=tzero-wire2Tzero/c2-pulse_width_range/4.0
  w_end=tzero-wire2Tzero/c2+pulse_width_range/4.0
  posw=where((wl ge w_start) and (wl le w_end),nposw)
  print,'nposw='+strtrim(string(nposw),2)
  ;; find locations of bins for mean base calculation
  b_start = out_of_pulse*time_step
  b_end = wl[nb-1]
  target_start = tzero + target_range/c2 - pulse_width_range
  target_end = tzero + target_range/c2 + pulse_width_range
  posb = where((wl ge b_start) and (wl le b_end) and ((wl le target_start) or (wl ge target_end)), nposb)

  ;; refine the window of target return such that it won't mess up wire
  ;; signal search of near ranges
  use_apriori_tzero = 0b
  if (target_start lt w_end) then begin
    print, 'target_start time = ', target_start
    print, 'w_end time = ', w_end
    target_start = tzero + target_range/c2 - pulse_width_range/4.0
    if (target_start lt w_end) then begin
      midtime = (target_start + w_end)/2.0
      print, 'midtime = ', midtime
      target_start = midtime
      ;; tell the program to use a priori tzero
      use_apriori_tzero = 1b
      ap_tzero = tzero
    endif else begin
      target_start = w_end
    endelse 
  endif 
  print, 'refined target_start time = ', target_start
  postarget = where((wl ge target_start) and (wl le target_end), npostarget)

  for i=0L,nl-1L do begin
    satmask=bytarr(ns)
    ;; it's a fake data cube from stationary scan, use all waveforms in
    ;; a scan line. 
    index=where((mask_all[*,i] ne 0), count)
    if (count gt 0L) then begin
      pointsz=long64(i)*long64(bufrs)
      data=read_binary(inlun,data_start=pointsz,data_dims=[ns,nb],data_type=type)
      maxtemp = max(data, dimension=2)
      sat_pos = where(maxtemp ge sat_test, count_sat)
      if (count_sat gt 0) then begin
        satmask[sat_pos]=1b
        index=where((mask_all[*,i] ne 0) and (satmask eq 0b), count)
        if (count le 0) then goto,nodata
      endif
      ;
      d = double((data[index,*]))
      data=0b
      n = long(n) + 1L
      ;; median waveform of this 'scan line'
      ;; due to severe noise in waveforms, median has been proved to be more
      ;; stable than simple mean. 
      temp = median(d, dimension=1, /even)
      ;; temp=total(d,1,/double)/double(count) 
      ;; median of waveform bins after the wire signal, it can give us the noise
      ;;base level.
      tmpmedian = median(d[*, posb], dimension=2, /even)
      ;; mean of medians of waveforms in this scan line as the mean base of this
      ;;scan line. 
      temp2 = mean(tmpmedian, /double)
      ;; temp2=total(d[*,out_of_pulse:nb-1],/double)/(double(count)*double(nb-out_of_pulse))
      ;
      nt=n_elements(temp)
      ;; get the max of wire return pulse
      tempmax=max(reform(temp[posw]),nct)
      mpos=posw[0]+nct
      line_scale[i]=target_dn/(tempmax-temp2)
      if (get_info_stats) then begin
        save[0,i]=float(i)
        save[1,i]=float(count)
        ;      save[2,i]=max(temp,mpos)-temp2
        save[2,i]=tempmax-temp2 ; raw DN max subtracts mean base DN
        save[3,i]=float(mpos)
        save[4,i]=float(temp2)  ; mean base DN from median of waveform bins
        save[5,i]=float(line_scale[i])
        save[6,i]=float(sqrt(total((temp-temp2)^2,/double)))
        ;      save[7,i]=(wl[1]-wl[0])*save[6,i]/(max(temp)-temp2)
        save[7,i]=(wl[1]-wl[0])*save[6,i]/save[2,i]
      endif
      ;; notice 'temp' here is from all waveforms and thus includes return
      ;;signals back from targets other than wire. Thus we first change the bins
      ;;in the expected target range to mean base
      temp[postarget] = temp2
      sum = sum + temp
      sum2 = sum2 + total(d^2, 1,/double)/double(count)
    endif else begin
      nodata:
      if (get_info_stats) then begin
        save[0,i]=float(i)
        save[1,i]=0.0
        save[2,i]=0.0
        save[3,i]=0.0
        save[4,i]=0.0
        save[5,i]=1.0
        save[6,i]=0.0
        save[7,i]=0.0
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
  
  if (n le 0) then begin
    print,'No valid casing values at all!!'
    err=66
    goto,cleanup
  endif
  d=double(n)
  ;; notice the pulse here includes return signals back from targets other than
  ;; wire. 
  pulse=float(sum/d)
  sig=float(sqrt(abs(sum2-sum^2/d)/d))
  ;help,sig
  ;make space again
  sum=0b
  sum2=0b
  d=0b
  
  ;mu and sig are length the number of bands
  ;; now the pulse and sig includes signals of all targets, not just the wire
  mean_base = median(pulse[posb], /double, /even)
  mean_base_sig = median(sig[posb], /double, /even)
  ;; mean_base=total(pulse[out_of_pulse:nb-1],/double)/double(nb-out_of_pulse)
  ;; mean_base_sig=total(sig[out_of_pulse:nb-1],/double)/double(nb-out_of_pulse)
  cv_base=100.0*mean_base_sig/mean_base
  
  pulse=pulse-mean_base
  ;casing_power=sqrt(total(pulse^2,/double))
  casing_power=total(pulse,/double)
  if (abs(casing_power) lt 1.0e-6) then casing_fwhm=0.0 else $
    casing_fwhm=(wl[1]-wl[0])*casing_power/max(pulse)
  ;
  ;=======================================================
  ;if get_info_stats set then write out the saved data
    
  o_path=file_dirname(out_satfix_name)
  if(get_info_stats) then begin
    ;strtrim(DWELCubeFile,2)
    n_base=strlen(out_satfix_name)
    n_dot=strpos(out_satfix_name,'.',/reverse_search)
    ;
    if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
      clog_file=strtrim(out_satfix_name,2)+'_nu_basefix_trace.log'
    endif else begin
      clog_file=strmid(strtrim(out_satfix_name,2),0,n_dot)+'_nu_basefix_trace.log'
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
      print,'DWEL_baseline_sat_fix_cmd terminating'
      print,' '
      err=7
      goto,cleanup
    endif
    ;
    printf,ctfile,strtrim('DWEL calibration casing trace Log File',2)
    printf,ctfile,strtrim('Run made at: '+systime(),2)
    printf,ctfile,strtrim('Input Cube File: '+strtrim(DWELCubeFile,2),2)
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
    ;; for stationary scan, do NOT scale so that we can trace back the real DN
    ;; for calibration estimate
    scale_mean = 1.0
    scale_cv = 0.0

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
      ;    pos=where(reform(mask_all[*,i]) ne 0,npos)
      ;    if (npos le 0) then begin
      ;      power=0.0
      ;    endif else begin
      ;      power=total(reform(anc_data[*,i,4]))/float(npos)
      ;    endelse
      outstring=strtrim(string(save[0,i],format='(f10.0)'),2)+','+ $
        strtrim(string(save[1,i],format='(f10.0)'),2)+','+ $
        strtrim(string(save[2,i],format='(f10.3)'),2)+','+ $
        strtrim(string(save[3,i],format='(f10.0)'),2)+','+ $
        strtrim(string(save[4,i],format='(f10.3)'),2)+','+ $
        strtrim(string(save[5,i],format='(f10.3)'),2)+','+ $
        strtrim(string(save[6,i],format='(f10.3)'),2)+','+ $
        strtrim(string(save[7,i],format='(f10.3)'),2)
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
  time = wl
  ;  p_time = time
  ;  pulse = sum / double(n)
  ;  sig = sqrt((sum2 / double(n) - pulse^2)*double(n)/double(n-1))
  CasingMeanWfMax = max(reform(pulse[posw]), nct)
  Tzero_I=posw[0]+nct
  print, 'Initial Tzero before baseline fix = ', time[Tzero_I], ' ns'
  tlow=time[Tzero_I] - 1.5*pulse_width_range
  thigh=time[Tzero_I] +1.5*pulse_width_range
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
  baseline[postarget] = 0.0
  
  pulse = pulse - baseline
  
  CasingMeanSig=sig[Tzero_I]
  CasingMeanCV=100.0*CasingMeanSig/CasingMeanWfMax

  ;; if (out_of_pulse-bins_toward_wire-Tzero_I)*time_step lt thigh then begin
  ;;   thigh = (out_of_pulse-bins_toward_wire-Tzero_I)*time_step
  ;; endif 
  tmpind = where((time ge tlow) and (time le thigh), count)
  if count gt 0 then begin
    casing_power=total(reform(pulse[tmpind]),/double)
    tmpmax = max(pulse[tmpind])
    if (abs(casing_power) lt 1.0e-6) then casing_fwhm=0.0 else $
      casing_fwhm=(wl[1]-wl[0])*casing_power/tmpmax
  endif else begin
    casing_power = 0.0
    casing_fwhm = 0.0
  endelse 
  tmpind=0b

  Tzero = time[Tzero_I]
  ;; interpolate peak location
  istat = peak_int(time[[Tzero_I-1, Tzero_I, Tzero_I+1]], pulse[[Tzero_I-1, Tzero_I, Tzero_I+1]], time_int, pulse_int, offset)
  Tzero = time_int
  print, 'Initial Tzero after interpolation = ', Tzero, ' ns'
  
  delta= wire2Tzero/c2
  print,'delta=',delta, ' ns'
  print,''
  
  Tzero=Tzero+delta
  print,'Shifted Tzero',Tzero, ' ns'

  if (use_apriori_tzero) then begin
    print, 'Target is too close and contaminates wire signal!'
    print, 'Now use the a priori Tzero (ns), ', ap_tzero
    tzero = ap_tzero
  endif 
  
  time=time-Tzero
  
  tend=(Tzero - pulse_width_range)>0.0
  tnice=(Tzero - pulse_width_range/2.0)>0.0
  shottz=fix(tnice/time_step)
  shotend=fix(tend/time_step)
  print,'tend=',tend
  print,'shotend,shottz=',shotend,shottz
  
  ;=================================================================================================
  ;write out the information into a spectral library to check later
  write_baseline_file:
  
  ;Determine library file name
  n_base=strlen(out_satfix_name)
  n_dot=strpos(out_satfix_name,'.',/reverse_search)
  if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
    baseline_library=strtrim(out_satfix_name,2)+'_base_Lib.sli'
  endif else begin
    baseline_library=strmid(out_satfix_name,0,n_dot)+'_base_Lib.sli'
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
  free_lun,lun
  
  ;set up the ENVI header for the output library
  descrip='Baseline Fit Information from: '+strtrim(DWELCubeFile,2)
  envi_setup_head,fname=baseline_library,ns=nb,nl=3,nb=1,$
    xstart=0,ystart=0,file_type=4,interleave=0,$
    data_type=sense.type,spec_names=lib_nam,$
    descrip=descrip, wl=time,$
    zplot_titles=['Time','Value'], $
    /write
    
  envi_open_file,baseline_library,r_fid=lib_fid,/no_interactive_query,/no_realize
  envi_file_mng,id=lib_fid,/remove
  
  out_mat=0b
  baseline_corr=0b
  ;=================================================================================================
  
  ;set up the DWEL header information for the base fixing
  ;
  
  dwel_pulse_model_dual_nsf, wavelength, i_val, t_val, r_val, p_range, p_time, pulse_model
  pulse_model_name = 'NSF_DWEL_Pulse_Model'
  ;; if the input data is from  Oz DWEL, keopsys lasers, 
  if strcmp(laser_man, 'keopsys', /fold_case) then begin
    DWEL_pulse_model_dual_oz, wavelength, i_val, t_val, r_val, p_range, p_time, pulse_model
    pulse_model_name = 'Oz_DWEL_Pulse_Model'
  endif 

  p_troughloc = i_val[3]
  p_scdpeakloc = i_val[4]
  
  model_fwhm=total(pulse_model)*time_step
  
  DWEL_base_fix_info=strarr(21)
  DWEL_base_fix_info=[$
    'Program=dwel_baseline_sat_fix_cmd_nsf',$
    'Descr=DWEL New Base Fix Settings with casing power',$
    'Processing Date Time='+strtrim(systime(),2),$
    'Pulse='+pulse_model_name,$
    'Comment=Tzero is the time at which the peak of the output pulse occurs',$
    'Tzero='+strtrim(string(Tzero,format='(f10.3)'),2),$
    'srate='+strtrim(string(srate,format='(f10.2)'),2),$
    'out_of_pulse='+strtrim(string(out_of_pulse,format='(i10)'),2),$
    'Target_DN='+strtrim(string(target_dn,format='(f10.2)'),2),$
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
  DWEL_base_fix_info=strtrim(DWEL_base_fix_info,2)
  
  ;=====================================================================
  ;now apply the correction to the image
  
  nb_out=nb
  nl_out=nl
  ns_out=ns
  wl_range=c2*time
  
  ;Open the output file for BIL tiling
  text_err=0
  
  ;=====================================================================
  openw, osatfile, out_satfix_name,/get_lun,error=text_err
  if (text_err ne 0) then begin
    print, strtrim('Halting DWEL_static_wire_baseline_sat_fix', 2)
    print, strtrim(['Error opening output file '+strtrim(out_satfix_name,2)], 2)
    err=8
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
  
  runin_mask=fltarr(nb_out)+1.0
  numpos=indgen(shottz-shotend+1)+shotend
  ord=float(numpos-shotend)/float(shottz-shotend)
  
  print,'numpos ends=',numpos[0],numpos[shottz-shotend]
  help,ord
  
  print,'numpos'
  print,numpos
  print,''
  print,'ord'
  print,ord
  
  merge=(ord^2)*(3.0-2.0*ord)
  print, 'merge'
  print, merge
  runin_mask[0:shotend-1]=0.0
  runin_mask[shotend:shottz]=merge
  runin_mask[shottz+1:nb_out-1]=1.0

  ;; ;; if a target range is provided, it is interpreted as the maximum possible
  ;; ;; range of the targe in the static scan. 
  ;; ;; suppress all waveform bins after this range to zeros to avoid complications
  ;; ;;in later processing caused by noise. 
  ;; if n_elements(target_range) ne 0 or arg_present(target_range) then begin
  ;;   target_time = Tzero + double(target_range)/c2
  ;;   target_end = fix((target_time + pulse_width_range)/time_step)
  ;;   target_halfend = fix((target_time + pulse_width_range/2.0)/time_step)
  ;;   runin_mask[target_end+1:nb_out-1] = 0.0
  ;;   numpos=indgen(target_end-target_halfend+1) + target_halfend
  ;;   ord=float(numpos-target_halfend)/float(target_end-target_halfend)
  ;;   merge = ((1.0-ord)^2)*(1.0+2.0*ord)
  ;;   runin_mask[target_halfend:target_end] = merge
  ;;   target_start = fix((target_time - pulse_width_range)/time_step)
  ;;   wire_end = fix((Tzero + pulse_width_range)/time_step)
  ;;   if wire_end+1 le target_start-1 then begin
  ;;     runin_mask[wire_end+1:target_start-1] = 0.0
  ;;   endif 
  ;; endif 

  ;  help,runin_mask
  ;  mean_image=fltarr(ns_out,nl_out)
  ;  max_image=fltarr(ns_out,nl_out)
  ;==========================================
  sat_mean_image=fltarr(ns_out,nl_out)
  sat_max_image=fltarr(ns_out,nl_out)
  satfixeddata = fltarr(ns_out,nb_out)
  ;==========================================
  pos_pos=where(wl_range gt 0.0,npos_r)
  if(npos_r le 0) then begin
    print,'Undefined time range, npos_r=',npos_r
    print,'there were NO ranges gt 0.0!'
    err=10
    goto,cleanup
  endif
  fln=1.0/float(npos_r)
  
  ;print,'npos_r=',npos_r
  ;help,pos_pos
  
  sat_mask=bytarr(ns_out,nl_out)
  sun_sensor=lonarr(ns_out,nl_out)
  
  alltotbad=0L
  totsat=0L
  
  ;set up the pointer to read the input file again
  bufrs=long(nbytes)*long(nb_out)*long(ns_out)
  pointsz=long64(0)
  
  ;do the processing over the lines for BIL structure
  
  for i=0, nl_out-1 do begin
    ;first get the data tile
    data=read_binary(inlun,data_start=pointsz,data_dims=[ns_out,nb_out],data_type=type)
    pointsz=long64(pointsz)+long64(bufrs)
    temp=float(data)
    ;============================================
    ;saturation detection
    ;check if the waveform maximum is equal or larger sat_test=1023. If so the waveform is identified saturated
    maxtemp = max(temp, dimension=2)
    sat_pos = where(maxtemp ge sat_test, count_sat)
    locbad=0L
    ;check for REALLY bad places where there is total saturation that persists into later parts of the scan
    if (count_sat gt 0) then begin
      for kk=0,count_sat-1 do begin
        pover=where(temp[sat_pos[kk],*] ge sat_test,num_over)
        if (num_over gt 20) then begin
          temp[sat_pos[kk],*]=0.0
          mask_all[sat_pos[kk],i]=0
          alltotbad=long(alltotbad)+1L
          locbad=long(locbad)+1L
        endif
        pover=0b
      endfor
    endif
    if (locbad gt 0L) then begin
      maxtemp = max(temp, dimension=2)
      sat_pos = where(maxtemp ge sat_test, count_sat)
    endif
    
    ;first remove the constant baseline
    temp=temp-transpose(cmreplicate(mean_base,nb))##one_ns
    pos_z=where(mask_all[*,i] eq 0,count_z)
    if (count_z gt 0L) then begin
      temp[pos_z,*]=0.0
    endif
    
    ;this has been done to make satfix correct as data input - ie constant baseline removed
    ;============================================
    satfixeddata=temp
    ;saturation fix, replace the saturated waveform part with pulse model and corresponding peak intensity
    if (count_sat gt 0) then begin
      satflag=1
      totsat=long(totsat)+long(count_sat)
      totbad=0L
      for si=0,count_sat-1 do begin
        unsattemp=reform(temp[sat_pos[si], *])
        ;; record the Tzero pulse from the wire
        ;; find the wire signal maximum
        wirewf = unsattemp
        satflag = apply_sat_fix_nsf(reform(temp[sat_pos[si], *]), pulse_model, i_val, scale_mean, satfixedwf=unsattemp)
        if (satflag lt 2) then sat_mask[sat_pos[si],i]=1
        if (satflag le 0) then begin
          totbad=long(totbad)+1L
          satfixeddata[sat_pos[si], *]=0.0
          mask_all[sat_pos[si],i]=0
          sat_mask[sat_pos[si],i]=2
          ;          print, 'line=', i, ', sample=', sat_pos[si]
          continue
        endif else begin
          ;; if the substitute section by sat fix overlap with wire pulse, then
          ;; paste the wire signal back
          tmpmax = max(unsattemp, tmpind)
          if (tmpind - i_val[2] + i_val[1] - 1) le out_of_pulse then begin
            wiremax = max(wirewf[posw], wiremaxind)
            wiremaxind = wiremaxind + posw[0]
            unsattemp[wiremaxind-i_val[2]+i_val[1]:wiremaxind-i_val[2]+i_val[5]] $
              = unsattemp[wiremaxind-i_val[2]+i_val[1]:wiremaxind-i_val[2]+i_val[5]] $
              + pulse_model[i_val[1]:i_val[5]]*wiremax
          endif 
          satfixeddata[sat_pos[si], *] = unsattemp
        endelse
        unsattemp=0b
      endfor
      alltotbad=long(alltotbad)+long(totbad)
    endif
    temp=satfixeddata
    satfixeddata=0b
    
    ;======================================================================
    temp=temp-transpose(baseline)##one_ns
    temp=scale_mean*temp*(transpose(runin_mask)##one_ns)
        
    ;mask may have changed now
    pos_z=where(mask_all[*,i] eq 0,count_z)
    if (count_z gt 0L) then begin
      temp[pos_z,*]=0.0
    endif
    
    ;now find the bad data such as baseline sun brightness
    tmat=reform(temp[*,far_range:nb-1])
    test=total(tmat,2)/float(nb-far_range)
    test2=total(tmat^2,2)/float(nb-far_range)
    test2=sqrt(test2-test^2)
    test=test/test2
    tmat=0b
    test2=0b
    pos_sun=where(abs(test) ge sun_thresh,npos_sun)
    if (npos_sun gt 0) then begin
      mask_all[pos_sun,i]=0b
      temp[pos_sun,*]=0.0
      sun_sensor[pos_sun,i]=1L
    endif
    
    ;    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;replaced code
    ;    ; check if the mean of the first 100 waveform bins is abnormal
    ;    tmpmean = mean(temp[*, 0:99], dimension=2)
    ;    abnormalpos = where(tmpmean gt 10, tmpcount)
    ;    if (tmpcount gt 0) then begin
    ;      temp[abnormalpos, *] = temp[abnormalpos, *] - cmreplicate(tmpmean[abnormalpos], [1, nb_out])
    ;    endif
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    ;The satfix section has been taken from here to a previous spot
    
    ;=================================================
    ;round if integer else set back in data
    if (type lt 4 or type gt 9) then begin
      temp=fix(round(temp), type=2)
    endif
    
    ;write out the resulting tile
    writeu,osatfile,temp
    sat_mean_image[*,i]=fln*total(temp[*,pos_pos],2)
    sat_max_image[*,i]=float(max(temp[*,pos_pos],DIMENSION=2))
    ;=================================================
    
    data=0
    satfixeddata=0b
    temp=0b
  endfor
  
  free_lun,osatfile,/force
  free_lun,inlun,/force
  
  print,'number saturated=',totsat
  print,'number of bad=',alltotbad
  ;set the output data type
  if (type lt 4 or type gt 9) then begin
    out_type=2
  endif else begin
    out_type=4
  endelse
  
  ;clear up and complete the action
  ;
  data=0b
  temp=0b
  one_ns=0b
  one_nb=0b
  ;  envi_file_mng,id=infile_fid,/remove
  ;===================================================
  ;get names
  out_base=file_basename(out_satfix_name)
  ;get output_ancillary file name
  ;Set up the ancillary file
  dot = strpos(out_satfix_name,'.',/reverse_search)
  if ((dot lt 0) or ((strlen(out_satfix_name)-dot-1) ne 3)) then dot = strlen(out_satfix_name)
  ancfile = strmid(out_satfix_name, 0, dot)+'_ancillary.img'
  last=strpos(ancfile,path_sep(),/reverse_search)
  anc_base=file_basename(ancfile)
  
  ;===================================================
  ; Write saturation fixing parameters to header
  DWEL_sat_info=strarr(10)
  DWEL_sat_info = [$
    'Program=dwel_baseline_sat_fix_cmd_nsf',$
    'Title=Parameters for Saturation Fixing',$
    'Model='+strtrim(pulse_model_name,2),$
    'Sat test value='+strtrim(sat_test,2), $
    'Saturated_(pixels)='+strtrim(totsat,2), $
    'Bad_(pixels)='+strtrim(alltotbad,2), $
    'Stats_Format=(Num,Min,Max,Mean,RMS)',$
    'Range_Stats=('+'0'+')',$
    'Sat Fixed File='+out_base, $
    'Updated Ancillary File='+strtrim(anc_base,2)]
  DWEL_sat_info=strtrim(DWEL_sat_info,2)
  
  DWEL_sat_info=[DWEL_sat_info,'Sat_Fixed_File='+strtrim(out_base,2)]
  
  ;write out header for the output file
  descrip='DWEL base fix and Sat Fix applied to '+strtrim(out_base,2)
  band_names=strarr(nb_out)+'Range_Sample_(m)_'
  band_names=band_names+strtrim(string(indgen(nb_out)+1),2)+'_satfix'
  
  envi_setup_head,fname=out_satfix_name,ns=ns_out,nl=nl_out,nb=nb_out,$
    xstart=xstart+dims[1],ystart=ystart+dims[3],$
    data_type=out_type, interleave=ft_out, $
    wl=wl_range, inherit=inherit, $
    bnames=band_names,descrip=descrip, $
    zplot_titles=['Range (m)','Intensity'], $
    /write
    
  envi_open_file,out_satfix_name,r_fid=out_fid,/no_interactive_query,/no_realize
  
  ;write out the previous header records
  status=DWEL_put_headers(out_fid,DWEL_headers)
  ;
  ;write the new header(s) into the HDR file
  envi_assign_header_value, fid=out_fid, keyword='DWEL_base_fix_info', $
    value=DWEL_base_fix_info
  envi_assign_header_value, fid=out_fid, keyword='DWEL_sat_info', $
    value=DWEL_sat_info
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
    print, strtrim('Halting DWEL_baseline_sat_fix', 2)
    print, strtrim(['Error opening output file '+strtrim(ancname,2)], 2)
    err=11
    goto, cleanup
  endif
  
  pos_pos=where(mask_all ne 0)
  sat_meanmean=mean(sat_mean_image[pos_pos])
  sat_stddevmean=stddev(sat_mean_image[pos_pos])
  sat_mlow=sat_meanmean-4.0*sat_stddevmean
  sat_mhigh=sat_meanmean+4.0*sat_stddevmean
  sat_minmean=min(sat_mean_image[pos_pos]) > sat_mlow
  sat_maxmean=max(sat_mean_image[pos_pos]) < sat_mhigh
  sat_mean_image[pos_pos]=((4095.0*(sat_mean_image[pos_pos]-sat_minmean)/(sat_maxmean-sat_minmean) > 0.0) < 4095.0)
  
  print,'sat fixed data scale range=',sat_minmean,sat_maxmean
  writeu,ofile,long(sat_mask)
  writeu,ofile,long(sun_sensor)
  for j=2,4 do begin
    writeu,ofile,long(anc_data[*,*,j])
  endfor
  writeu,ofile,round(sat_max_image)
  ;  writeu,ofile,round(sat_mean_image)
  writeu,ofile,long(mask_all)
  writeu,ofile,round(100.0*zeniths)
  writeu,ofile,round(100.0*azimuths)
  free_lun, ofile,/force
  anc_data=0b
  sat_mean_image=0b
  
  envi_setup_head, fname=ancfile, $
    ns=ns_out, nl=nl_out, nb=9, $
    interleave=0, data_type=3, $
    /write, $
    bnames=['Sat Mask','Sun Mask','Scan Encoder','Rotary Encoder', $
    'Laser Power','Waveform Mean','Mask','Zenith','Azimuth']
    
  envi_open_file,ancfile,r_fid=anc_fid,/no_interactive_query,/no_realize
  
  ;write out the previous header records
  status=DWEL_put_headers(anc_fid,DWEL_headers)
  envi_assign_header_value, fid=anc_fid, keyword='DWEL_Adaptation', $
    value=DWEL_Adaptation
  ;; write new header records
  envi_assign_header_value, fid=anc_fid, keyword='DWEL_base_fix_info', $
    value=DWEL_base_fix_info
  envi_assign_header_value, fid=anc_fid, keyword='DWEL_sat_info', $
    value=DWEL_sat_info
  envi_write_file_header, anc_fid
  
  envi_file_mng,id=anc_fid,/remove
  ;=====================================================================
  ;
  sat_mask=0b
  ;=====================================================================
  
  ; Do the final cleanup
  cleanup:
  free_lun,lun,/force
  free_lun,inlun,/force
  free_lun,ofile,/force
  free_lun,tfile,/force
  free_lun,ctfile,/force
  free_lun,osatfile,/force
  if (err ne 0) then print,'Returning from DWEL_Baseline_Sat_Fix_Cmd_nsf with error'
  heap_gc,/verbose
  ;

  ;; write processing time summary
  print, '**************************************************'
  print, 'Processing program = dwel_static_wire_baseline_sat_fix_cmd_nsf'
  print, 'Input DWEL cube file size = ' + $
    strtrim(string(double(procfilesize)/(1024.0*1024.0*1024.0)), 2) + ' G'
  print, 'Processing time = ' + strtrim(string((systime(1) - starttime)), 2) + ' ' + $
    'seconds'
  print, '**************************************************'

  print, '************************************************************'
  print, 'The minimum range of stationary scans this processing, '
  print, 'dwel_static_wire_baseline_sat_fix_cmd_nsf can accept = 5.0 m'
  print, 'If your stationary scan is from a closer target, '
  print, 'the processing may give unexpected results'
  print, '************************************************************'

  return
end
