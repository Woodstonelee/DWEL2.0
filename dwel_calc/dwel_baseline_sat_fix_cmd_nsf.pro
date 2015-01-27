;; baseline fix and saturation fix of NSF DWEL data
;; Zhan Li, zhanli86@bu.edu
;; Created in 2013 by Zhan Li
;; Last modified: 20140603 by Zhan Li

pro DWEL_Baseline_Sat_Fix_cmd_nsf, DWELCubeFile, ancillaryfile_name, out_satfix_name, Casing_Range, Casing_Type, get_info_stats, zen_tweak, err, wire=wire, settings=settings
;+
;PURPOSE:
;; Fix baseline and saturation of a DWEL data cube. 
;; 
;INPUTS:
;; DWELCubeFile = string, the full file name of the DWEL cube file
;; 
;; ancillaryfile_name = string, the full file name of the ancillary file of the DWEL
;; cube file. 
;; 
;; out_satfix_name = string, the full file name of the output baseline fixed and
;; saturation fixed waveform data. 
;; 
;; Casing_Range = two-element vector of numericals, [min_zen_angle,
;; max_zen_angle], the zenith range to search casing returns for inital Tzero,
;; laser power variation, and background noise level. 
;; 
;; Casing_Type = string, the type of the casing area given by
;; Casing_Range. Available input strings are, 
;; 'LAM', the on-board lambertian target. 
;; 'CASE', the edge of the instrument case. It is used usually when the
;; lambertian target is not on board or the cap of the lambertian target is not
;; removed during scanning. 
;; 'CAP', the cap of the lambertian target. This is an optional choice when the
;; cap of the lambertian target is not removed during scanning. 
;;
;; get_info_stats = byte integer, tell the program whether to collect and output
;; detailed stats of baseline and saturation fix. 
;;
;; zen_tweak = numerical, tweak the encoder of zenith point to refine the
;; calculation of zenith angle from scan encoder values. 
;; 
;OUTPUTS:
;; err = byte integer, return the error code of program running. 
;;
;KEYWORDS:
;; wire = binary keyword, tell the program whether the input data is collected
;; with wire or without wire. If set, with wire. If false, without wire. 
;;
;; settings = structure, one pair of tag and value gives a user-defined value of
;; one setting. The available tag names of settings and their values are, 
;; 'out_of_pulse', integer, the bin location where you are free of casing
;; returns in a waveform. It is used to estimate mean baseline. Default value =
;; 200. 
;; THE FOLLOWING TWO SETTINGS ARE REMOVED AND NOT AVAILABLE ANY MORE.
;; 'casing2Tzero', numerical, the distance between the scan mirror center and
;; the casing area provided by Casing_Range. If your Casing_Range is NOT
;; centered at nadir (upper limit is not 180 degrees), you MUST update this
;; setting with a correct value. Default, no setting.  
;; 'wire2Tzero', numerical, the distance between the scan mirror center and the
;; wire if the wire is present. This value is not actually used in this step but
;; will be recorded in the header file and starts to be used in the next step of
;; filtered post-fix processing. Default, no setting.
;;
;RETURNS:
;; None.
;;
;MODIFICATION HISTORY:
;; Created, 201406, by Zhan Li, zhanli86@bu.edu.
;; Changed and contributed heavily, 201411 - 20141219, by David Jupp. 
;-
;

  compile_opt idl2

  resolve_routine, 'DWEL_GET_HEADERS', /compile_full_file, /either
  resolve_routine, 'DWEL_SET_THETA_PHI_NSF', /compile_full_file, /either
  resolve_routine, 'DWEL_PULSE_MODEL_DUAL_NSF', /compile_full_file, /either
  resolve_routine, 'DWEL_PUT_HEADERS', /compile_full_file, /either
  resolve_routine, 'APPLY_SAT_FIX_NSF', /compile_full_file, /either
  resolve_routine, 'CMREPLICATE', /compile_full_file, /either
  resolve_routine, 'CLEAN_LINE', /compile_full_file, /either

  ;; get the size of input file to be processed. It will be used in later
  ;; summary of processing time. 
  procfilesize = file_info(DWELCubeFile)
  procfilesize = procfilesize.size
  ;; get the time now as the start of processing
  starttime = systime(1)

  ;
  lun=99
  wlun=107
  inlun=105
  ofile=101
  osatfile=102
  tfile=98
  ctfile=35
  fname=''
  o_name=''
  ctfile=30
  err=0
  
  ;======================================================================
  ;more internal settings
  ;always get the pulse information and write out files at this time
  get_info_stats=1b
  ;skip these values - currently noisy and useless, in unit of bins
  before_casing=100
  ;target dn preset but set later
  target_dn = 512.0
  if strcmp(Casing_Type, 'CAP', /fold_case) then begin
    target_dn = 150.0
  endif
  if strcmp(Casing_Type, 'CASE', /fold_case) then begin
    target_dn = 512.0 ;; WRONG!!! need update
  endif 
  ;used in sun background check to be away from significant return pulses
  far_range=650
  ;threshold on raise of baseline by sun background radiation, note these are
  ;scaled and initially basefixed units 
  sun_thresh=1.5
  ;the "FWHM" of outgoing pulse (now just a number used and calibrated)
  outgoing_fwhm = 5.1
  ;the full width of outgoing pulse where intensity is below 0.01 of maximum, in
  ;unit of ns. 
  pulse_width_range = 40.0
  ;saturation test value in DN
  sat_test=1023L
  ;set the 5 point gaussian filter
  filtg5=set_gfilt(2)

  finalsettings = { bsfixsettings, $
    ; position where you are free of the casing effects - for mean baseline
    out_of_pulse:200}
    ;; ;distance from casing (edge of casing) to the true Tzero position at mirror,
    ;; ;unit=metres 
    ;; casing2Tzero:0.055}
    ;; ; distance between scan mirror center and the wire if present, unit=meters. 
    ;; wire2Tzero:0.299}
  ;; if user provides settings
  if n_elements(settings) ne 0 or arg_present(settings) then begin
    finalsettings = update_struct_settings(finalsettings, settings)
  endif 
  out_of_pulse = finalsettings.out_of_pulse
  ;; casing2Tzero = finalsettings.casing2Tzero
  ;; wire2Tzero = finalsettings.wire2Tzero

  casing2Tzero = 0.055 ; unit, meter, default is for NSF

  ;======================================================================
  
  print,'pulse_width_range='+strtrim(string(pulse_width_range),2)
  
  ;clean up any fids which are no longer where they were!
  ;ENVI issue that is annoying and leads to confusion
  clean_envi_file_fids
  
  print,'entering basefix satfix program'
  
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
  n_base=strlen(fname)
  n_dot=strpos(fname,'.',/reverse_search)
  ;  if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
  ;    anc_name=strtrim(fname,2)+'_ancillary'
  ;  endif else begin
  ;    ancillaryfile_name=strmid(fname,0,n_dot)+'_ancillary.img'
  ;  endelse
  
  if(not file_test(ancillaryfile_name)) then begin
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
  status=DWEL_get_headers(infile_fid,DWEL_headers)
  
  if (not status) then begin
    print,strtrim('Bad FID in DWEL_get_headers! DWEL Header setup cancelled!',2)
    print,'Input File: '+strtrim(DWELCubeFile,2)
    err=5
    goto, cleanup
  endif
  
  ;help,dwel_headers,/structures
  
  if (DWEL_headers.headers_present le 0s or not DWEL_headers.run_present) then begin
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
  low=casing_range[0]
  high=casing_range[1]
  ;set wire limits
  w_zen_min=25.0
  w_zen_max=65.0
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
  if strcmp(laser_man, 'keopsys', /fold_case) then begin
    casing2Tzero = 0.065 ; unit, meter
  endif 
  
  ;input some planes of data from the ancillary file
  anc_data=lonarr(ns,nl,9)
  for j=0,8 do begin
    anc_data[*,*,j]=long(ENVI_GET_DATA(fid=ancillaryfile_fid, dims=[-1L,0,ns-1,0,nl-1], pos=j))
  endfor
  
  ;; scan encoder of zenith point, try to find it first from header file. If not
  ;; found, then the default values in dwel_set_phi_theta will be used later. 
  zenenc = -1
  ;; get some information for header file
  DWEL_Adaptation = ENVI_GET_HEADER_VALUE(ancillaryfile_fid, 'DWEL_Adaptation', undefined=undef)
  if undef then begin
    DWEL_Adaptation = ''
    if (strpos(DWEL_headers.f_base, '1064') ne -1) then begin
      wavelength = 1064
    endif
    if (strpos(DWEL_headers.f_base, '1548') ne -1) then begin
      wavelength = 1548
    endif
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

    match = -1
    for i=0,n_elements(info)-1 do begin
      if (strmatch(info[i],'*Scan encoder of zenith point*', /fold_case)) then match=i
    endfor
    if match ge 0 then begin
      text=strtrim(info[match],2)
      print,'text=',text
      k=strpos(text,'=')
      print,'extract=',strtrim(strmid(text,k+1),2)
      zenenc=fix(strtrim(strmid(text,k+1),2), type=3)
    endif
    print, 'Extracted ZenEnc = '+strtrim(string(zenenc), 2)
  endelse
  
  if (wavelength eq 1064) then begin
    target_dn = 512.0
    if strcmp(Casing_Type, 'CAP', /fold_case) then begin
      target_dn = 150.0
    endif
    if strcmp(Casing_Type, 'CASE', /fold_case) then begin
      target_dn = 512.0 ;; WRONG!!! need update
    endif 
  endif else begin
    target_dn = 509.0
    if strcmp(Casing_Type, 'CAP', /fold_case) then begin
      target_dn = 75.0
    endif
    if strcmp(Casing_Type, 'CASE', /fold_case) then begin
      target_dn = 509.0 ;; WRONG!!! need update
    endif 
  endelse
  
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
  
  if zenenc eq -1 then begin
    ;set up a structure and push it onto the heap
    sav={ $
      Nshots:ns,$
      Nscans:nl,$
      ShotZen:scan_encoder,$
      ShotAzim:rotary_encoder $
      }
  endif else begin
    ;set up a structure and push it onto the heap
    sav={ $
      Nshots:ns,$
      Nscans:nl,$
      ShotZen:scan_encoder,$
      ShotAzim:rotary_encoder, $
      ZenEnc:zenenc $
      }    
  endelse 
    
  ;now put the data on the heap with a pointer
  p_stat=ptr_new(sav,/no_copy)
  ;compute the zeniths and azimuths
  
  status = dwel_set_theta_phi_nsf(p_stat,zen_tweak)
  ;put the results into the local arrays
  zeniths=(*p_stat).ShotZen
  azimuths=(*p_stat).ShotAzim
  ptr_free, p_stat
  
  ;------------------------------------------------------------------------------
  ;now get the output file name
  
  output:
  
  n_base=strlen(fname)
  n_dot=strpos(fname,'.',/reverse_search)
  
  ;======================================================================
  
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
  
  ;set up the pulse model
  ;; default pulse model is from NSF DWEL, manlight lasers.
  DWEL_pulse_model_dual_nsf, wavelength, i_val, t_val, r_val, p_range, p_time, pulse_model
  pulse_model_name = 'NSF_DWEL_Pulse_Model'
  ;; if the input data is from  Oz DWEL, keopsys lasers, 
  if strcmp(laser_man, 'keopsys', /fold_case) then begin
    DWEL_pulse_model_dual_oz, wavelength, i_val, t_val, r_val, p_range, p_time, pulse_model
    pulse_model_name = 'Oz_DWEL_Pulse_Model'
  endif 

  model_fwhm=total(pulse_model)*time_step
  
  p_index=p_time/time_step
  help,p_index
  ;print,'p_index='+strjoin(strtrim(p_index,2),',',/single)+']'
  ;print,'i_val=['+strjoin(strtrim(i_val,2),',',/single)+']'
  
  ;get mean pulse and baseline from the given casing area designated
  ;by the zenith angles.
  n = long(0)
  sum = dblarr(nb)
  sum2 = dblarr(nb)
  temp=dblarr(nb)
  line_scale=fltarr(nl)+1.0
  line_offset=dblarr(nl)
  casing_max=fltarr(nl)
  casing_num=lonarr(nl)
  if (get_info_stats) then save=fltarr(8,nl)
  count=long(0)
  ;set up the pointer for read_binary
  bufrs=long(nbytes)*long(nb)*long(ns)
  pointsz=long64(0)
  ;
  print,'bufrs='+strtrim(string(bufrs),2)
  casing_bad_int=0L
  ;
  for i=0L,nl-1L do begin
    satmask=bytarr(ns)
    index=where((mask_all[*,i] gt 0) and ((zeniths[*,i] ge Casing_Range[0]) and (zeniths[*,i] le Casing_Range[1])), count)
    if (count lt 50L) then begin
      count=0L
      index=0b
      goto,nodata
    endif
    if (count gt 0L) then begin
      ;      data = envi_get_slice(fid=infile_fid, line=i, /bil)
      pointsz=long64(i)*long64(bufrs)
      data=read_binary(inlun,data_start=pointsz,data_dims=[ns,nb],data_type=type)
      maxtemp = max(data, dimension=2)
      sat_pos = where(maxtemp ge sat_test, count_sat)
      if (count_sat gt 0) then begin
        satmask[sat_pos]=1b
        index=where((mask_all[*,i] ne 0) and ((zeniths[*,i] ge Casing_Range[0]) and (zeniths[*,i] le Casing_Range[1]) $
          and (satmask eq 0b)), count)
        if (count le 50L) then goto,nodata
      endif
      ;
      d = double(reform(data[index,*]))
      data=0b
      n = long(n) + 1L
      temp=total(d,1,/double)/double(count)
      temp2=total(d[*,out_of_pulse:nb-1],/double)/(double(count)*double(nb-out_of_pulse))
      ;
      nt=n_elements(temp)
      tempmax=max(reform(temp[before_casing:nt-1]),nct)
      mpos=before_casing+nct
      ; interpolate peak location
      istat = peak_int([float(mpos-1), float(mpos), float(mpos+1)], temp[[mpos-1, mpos, mpos+1]], truepos, value, offset)
      if (istat ne 0) then begin 
        casing_bad_int=casing_bad_int+1L
        casing_max[i]=tempmax-temp2
      endif else begin
        casing_max[i] = value - temp2
      endelse 
      casing_max[i] = tempmax-temp2
      ;  print,'truepos='+strtrim(string(truepos),2)
      casing_num[i]=count      
      line_scale[i]=target_dn/(tempmax-temp2)
      line_offset[i]=temp2
      if (get_info_stats) then begin
        save[0,i]=float(i)
        save[1,i]=float(count)
        ;      save[2,i]=max(temp,mpos)-temp2
        save[2,i]=tempmax-temp2
        save[3,i]=float(truepos)
        save[4,i]=float(temp2)
        save[5,i]=float(line_scale[i])
        save[6,i]=float((total((temp-temp2)[mpos-80:mpos+80],/double)))
        ;      save[7,i]=(wl[1]-wl[0])*save[6,i]/(max(temp)-temp2)
        save[7,i]=(wl[1]-wl[0])*save[6,i]/save[2,i]
      endif
      sum = sum + temp
      sum2 = sum2 + total(d^2, 1,/double)/double(count)
    endif else begin
      nodata:
      count=0L
      index=0b
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
  
  if (casing_bad_int gt 0L) then print,'There were '+strtrim(string(bad_int),2)+' bad casing interpolations! ****'
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
  pulse=float(sum/d)
  sig=float(sqrt(abs(sum2-sum^2/d)/d))
  ;help,sig
  ;make space again
  sum=0b
  sum2=0b
  d=0b
  
  ;mu and sig are length the number of bands
  mean_base=total(pulse[out_of_pulse:nb-1],/double)/double(nb-out_of_pulse)
  mean_base_sig=total(sig[out_of_pulse:nb-1],/double)/double(nb-out_of_pulse)
  cv_base=100.0*mean_base_sig/mean_base
  
  pulse=pulse-mean_base
  
  nt=n_elements(pulse)
  maxpulse=max(reform(pulse[before_casing:nt-1]),nct)
  kmax=before_casing+nct
    
  help,pulse
  ;maxpulse=max(pulse,kmax)
  print,'kmax='+strtrim(string(kmax),2)
  print,'kmax+80='+strtrim(string(kmax+80),2)
  print,'kmax-80='+strtrim(string(kmax-80),2)
  
  ;casing_power=sqrt(total(pulse^2,/double))
  casing_power=total(pulse[kmax-80:kmax+80],/double)
  if (abs(casing_power) lt 1.0e-6) then casing_fwhm=0.0 else $
    casing_fwhm=(wl[1]-wl[0])*casing_power/maxpulse
  ;
  ;set up line offset so that offset is removed by line for improved accuracy
  valid_line=bytarr(nl)+1b
  pos_base=where(casing_num le 0,npos_base)
  if (npos_base gt 0) then begin
    valid_line[pos_base]=0b
  endif
  sm_line_offset=fltarr(nl)
  clean_line,line_offset,valid_line,sm_line_offset,filtg5,err
  
  pos_base=where(valid_line,npos_base)
  mean_line_offset=mean(sm_line_offset[pos_base],/double)

  ;; ***************************************************************************
  ;; calculate a scaling factor for each scan line to correct laser power
  ;; variation. not using line_scale is because it is noisy. 
  ;; first smooth the return intensities from standard target over scan lines
  ;; with clean_line function which does gap filling, median filtering and
  ;; Gaussian filtering to remove any possible noise or deviant in the mean
  ;; target returns from lines. 
  sm_casing_max = fltarr(nl)
  clean_line, casing_max, valid_line, sm_casing_max, filtg5, err
  ;; get scaling factor line by line to correct laser power variation
  sm_line_scale = target_dn / sm_casing_max

  valid_line=0b
  ;======================================================================
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
      scale_mean=mean(sm_line_scale[pos_sav],/double)
      scale_cv=stddev(sm_line_scale[pos_sav],/double)
      scale_cv=100.0*scale_cv/scale_mean
      old_scale_mean=mean(line_scale[pos_sav],/double)
      old_scale_cv=stddev(line_scale[pos_sav],/double)
      old_scale_cv=100.0*old_scale_cv/old_scale_mean
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
      old_scale_mean=1.0
      old_scale_cv=0.0
    endelse
    printf,ctfile,'Stats,Mean,CV(%)'
    outstring='Base_Stats(DN)='+strtrim(string(mean_base,format='(f14.3)'),2)+','+ $
      strtrim(string(cv_base,format='(f14.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    outstring='Casing_Stats(DN)='+strtrim(string(casing_mean,format='(f14.3)'),2)+','+ $
      strtrim(string(casing_cv,format='(f14.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    outstring='Pos_Stats(Bin)='+strtrim(string(pos_mean,format='(f14.3)'),2)+','+ $
      strtrim(string(pos_cv,format='(f14.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    outstring='Offset_Stats(DN)='+strtrim(string(off_mean,format='(f14.3)'),2)+','+ $
      strtrim(string(off_cv,format='(f14.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    outstring='Scale_Stats='+strtrim(string(old_scale_mean,format='(f14.3)'),2)+','+ $
      strtrim(string(old_scale_cv,format='(f14.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    outstring='sm_Scale_Stats='+strtrim(string(scale_mean,format='(f14.3)'),2)+','+ $
      strtrim(string(scale_cv,format='(f14.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    outstring='Mean Casing_Total(DN)='+strtrim(string(casing_power,format='(f14.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    outstring='Casing_fwhm(ns)='+strtrim(string(casing_fwhm,format='(f14.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    printf,ctfile,'Line_Num,Casing_Num,Casing_Max_Val(DN),Casing_Max_Pos(Bin),Offset(DN),Line_Scale,Casing_Total(DN),Casing_FWHM(ns),sm_Casing_Max_Val(DN),sm_Line_Scale'
    flush, ctfile
    ;
    for i=0L,nl-1L do begin
      ;    pos=where(reform(mask_all[*,i]) ne 0,npos)
      ;    if (npos le 0) then begin
      ;      power=0.0
      ;    endif else begin
      ;      power=total(reform(anc_data[*,i,4]))/float(npos)
      ;    endelse
      outstring=strtrim(string(save[0,i],format='(f14.0)'),2)+','+ $
        strtrim(string(save[1,i],format='(f14.0)'),2)+','+ $
        strtrim(string(save[2,i],format='(f14.3)'),2)+','+ $
        strtrim(string(save[3,i],format='(f14.3)'),2)+','+ $
        strtrim(string(save[4,i],format='(f14.3)'),2)+','+ $
        strtrim(string(save[5,i],format='(f14.3)'),2)+','+ $
        strtrim(string(save[6,i],format='(f14.3)'),2)+','+ $
        strtrim(string(save[7,i],format='(f14.3)'),2)+','+ $
        strtrim(string(sm_casing_max[i],format='(f14.3)'),2)+','+ $
        strtrim(string(sm_line_scale[i],format='(f14.3)'),2)
      ;
      printf,ctfile,outstring
    endfor
    ;
    flush, ctfile
    free_lun, ctfile,/force
    save=0b
  endif
  ;======================================================================
  
  ; initial time from current data cube before baseline fix
  time = wl
  ;  p_time = time
  ;  pulse = sum / double(n)
  ;  sig = sqrt((sum2 / double(n) - pulse^2)*double(n)/double(n-1))
  CasingMaxWfMean=max(reform(pulse[before_casing:n_elements(pulse)-1]), nct)
  posnl=where(casing_num gt 0L,nposnl)
  CasingMeanWfMax=mean(casing_max[posnl],/double)
  casingMaxratio=CasingMeanWfMax/CasingMaxWfMean
  Tzero_I=before_casing+nct
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
  
  pulse = pulse - baseline
  
  ;  CasingMeanWfMax = max(pulse, Tzero_I)
  CasingMeanSig=sig[Tzero_I]
  CasingMeanCV=100.0*CasingMeanSig/CasingMeanWfMax
  
  ;casing_power=sqrt(total(pulse^2,/double))
  casing_power=total(reform(pulse[tmpind]),/double)
  if (abs(casing_power) lt 1.0e-6) then casing_fwhm=0.0 else $
    casing_fwhm=(wl[1]-wl[0])*casing_power/max(pulse)
  tmpind=0b
  Tzero=time[Tzero_I]
  print,'Initial Tzero after baseline fix = ',Tzero,' ns'
  
  ;; interpolate peak location
  istat = peak_int(time[[Tzero_I-1, Tzero_I, Tzero_I+1]], pulse[[Tzero_I-1, Tzero_I, Tzero_I+1]], time_int, pulse_int, offset)
  Tzero = time_int
  print, 'Initial Tzero after interpolation = ', Tzero, ' ns'
  
  delta= -casing2Tzero/c2 
  print,'delta=',delta, ' ns'
  print,''
  
  Tzero=Tzero+delta
  print,'Shifted Tzero (ns)='+strtrim(string(Tzero),2)
  
  time=time-Tzero
  
  tend=(Tzero- pulse_width_range)>0.0
  tnice=(Tzero- pulse_width_range/2.0)>0.0
  shottz=fix(tnice/time_step)
  shotend=fix(tend/time_step)
  print,'tend=',tend
  print,'shotend,shottz=',shotend,shottz
  
  ;now set the scale factor
  
  scale_mean=target_dn/CasingMeanWfMax
  
  ;======================================================================
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
  
  ;======================================================================
  if keyword_set(wire) then begin
    ;get mean wire pulse from the given casing area designated
    ;by the zenith angles and the estimated area of the shots.
    w_start=tzero-0.4-pulse_width_range/4.0
    w_end=tzero-0.4+pulse_width_range/4.0
    posw=where((wl ge w_start) and (wl le w_end),nposw)
    print,'nposw='+strtrim(string(nposw),2)
    n = long(0)
    sum = dblarr(nposw)
    sum2 = dblarr(nposw)
    temp=dblarr(nposw)
    if (get_info_stats) then save_w=fltarr(8,nl)
    count=long(0)
    ;set up the pointer for read_binary
    bufrs=long(nbytes)*long(nb)*long(ns)
    pointsz=long64(0)
    save_tpos=fltarr(nl)
    save_tmax=fltarr(nl)
    save_valid=bytarr(nl)
    ;
    for i=0L,nl-1L do begin
      satmask=bytarr(ns)
      index=where((mask_all[*,i] gt 0) and ((zeniths[*,i] ge w_zen_min) and (zeniths[*,i] le w_zen_max)), count)
      if (count lt 50L) then begin
        count=0L
        index=0b
        goto,wire_nodata
      endif
      if (count gt 0L) then begin
        ;
        pointsz=long64(i)*long64(bufrs)
        data=read_binary(inlun,data_start=pointsz,data_dims=[ns,nb],data_type=type)
        ;      data=float(reform(data[*,posw]))
        maxtemp = max(data, dimension=2)
        sat_pos = where(maxtemp ge sat_test, count_sat)
        if (count_sat gt 0) then begin
          satmask[sat_pos]=1b
          index=where((mask_all[*,i] gt 0) and ((zeniths[*,i] ge w_zen_min) and (zeniths[*,i] le w_zen_max) $
            and (satmask eq 0b)), count)
          if (count le 50L) then goto,wire_nodata
        endif
        ;
        d = double(data[index,*])
        d=reform(d[*,posw])
        data=0b
        n = long(n) + 1L
        temp=total(d,1,/double)/double(count)
        ;
        nt=n_elements(temp)
        tempmax=max(temp,nct)
        ;print,'first tempmax,nct=',tempmax,nct
        ; new code to get peak nearest centre
        r=(temp-sm_line_offset[i])/(tempmax-sm_line_offset[i])
        test_pos=[replicate(0.0,nposw),r,replicate(0.0,nposw)]
        test_pos=convol(test_pos,filtg5)
        dc=deriv(test_pos)
        test_pos=0b
        d2c=deriv(dc)
        dr=dc[nposw:2*nposw-1]
        d2r=d2c[nposw:2*nposw-1]
        dc=0b
        d2c=0b
        ; Check neighbourhoods of derivative and second derivative of correlation for peaks
        bs1 = shift(dr,1)
        bs2=shift(dr,2)
        bsm1 = shift(dr,-1)
        bsm2=shift(dr,-2)
        test1 = (bs1 gt 0.0001 and bs2 gt 0.0001)
        testm1 = (bsm1 lt -0.0001 and bsm2 lt -0.0001)
        test = (d2r lt -0.01) and (r gt 0.25)
        peaks = where(((test) and (test1 and testm1)),nump)
        bs1=0b
        bs2=0b
        bsm1=0b
        bsm2=0b
        obj=abs(float(peaks)-float((nposw+1)/2))
        postz=min(obj,npt)
        nct=peaks[npt]
        tempmax=temp[nct]
        ;now see if there are two adjacent likely max positions (bins) and
        ;select max one 
        if (nump gt 1) then begin
          if (npt eq 0) then begin
            if (peaks[npt+1] eq peaks[npt]+1) then begin
              if (temp[peaks[npt+1]] gt tempmax) then npt=npt+1
            endif
          endif else if (npt eq nump-1) then begin
            if (peaks[npt-1] eq peaks[npt]-1) then begin
              if (temp[peaks[npt-1]] gt tempmax) then npt=npt-1
            endif
          endif else begin
            if (peaks[npt+1] eq peaks[npt]+1) then begin
              if (temp[peaks[npt+1]] gt tempmax) then npt=npt+1
            endif else if (peaks[npt-1] eq peaks[npt]-1) then begin
              if (temp[peaks[npt-1]] gt tempmax) then npt=npt-1
            endif
          endelse
        endif
        nct=peaks[npt]
        tempmax=temp[nct]
        ;print,'output tempmax,nct=',tempmax,nct
        r=0b
        obj=0b
        truepos=float(nct)
        ; interpolate peak location
        if ((nct gt 0)and (nct lt nt-1)) then begin
          istat = peak_int([float(nct-1), float(nct), float(nct+1)], temp[[nct-1, nct, nct+1]], truepos, value, offset)
          if (istat ne 0) then begin
            print,'bad interpolation for wire pos at line '+strtrim(string(i),2)
            print,'nct='+strtrim(string(nct),2)
            print,'peaks=['+strjoin(strtrim(string(peaks),2),',',/single)+']'
          endif
          tempmax=(1.0d0*temp[nct-1]+2.0d0*temp[nct]+1.0d0*temp[nct+1])/4.0d0
        endif
        ;
        mpos=posw[0]+nct
        truepos=truepos+float(posw[0])
        save_tpos[i]=float(truepos)
        save_tmax[i]=float(tempmax-sm_line_offset[i])
        save_valid[i]=1b
        
        if (get_info_stats) then begin
          save_w[0,i]=float(i)
          save_w[1,i]=float(count)
          save_w[2,i]=save_tmax[i]
          save_w[3,i]=save_tpos[i]
          save_w[4,i]=float(sm_line_offset[i])
          save_w[5,i]=save_tmax[i]
          save_w[6,i]=float((total((temp-sm_line_offset[i]),/double)))
          save_w[7,i]=(wl[1]-wl[0])*save_w[6,i]/save_w[2,i]
        endif
        sum = sum + temp
        sum2 = sum2 + total(d^2, 1,/double)/double(count)
      endif else begin
        wire_nodata:
        count=0L
        index=0b
        if (get_info_stats) then begin
          save_w[0,i]=float(i)
          save_w[1,i]=0.0
          save_w[2,i]=0.0
          save_w[3,i]=0.0
          save_w[4,i]=0.0
          save_w[5,i]=1.0
          save_w[6,i]=0.0
          save_w[7,i]=0.0
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
      err=67
      goto,cleanup
    endif
    dw=double(n)
    wpulse=float(sum/dw)-mean_base
    wsig=float(sqrt(abs(sum2-sum^2/dw)/dw))
    
    maxmeanwpulse=max(wpulse,kmax)
    
    posmaxw=where(save_valid gt 0b,nposmaxw)
    meanmaxwpulse=mean(save_tmax[posmaxw],/double)
    meantwzero=mean(save_tpos[posmaxw],/double)*time_step
    ;
    wire_power=total(wpulse,/double)*time_step
    if (abs(wire_power) lt 1.0e-6) then wire_fwhm=0.0 else $
      wire_fwhm=(wl[1]-wl[0])*wire_power/maxmeanwpulse
      
    print,'wire_fwhm='+strtrim(wire_fwhm,2)
    
    ;; get smoothed wire max and see how it looks like. 
    valid_line=bytarr(nl)+1b
    pos_base=where(save_valid le 0,npos_base)
    if (npos_base gt 0) then begin
      valid_line[pos_base]=0b
    endif
    sm_wire_max = fltarr(nl)
    clean_line, save_tmax, valid_line, sm_wire_max, filtg5, err

    ;make space again
    sum=0b
    sum2=0b
    d=0b
    
    ;if get_info_stats set then write out the saved data to stats file
    ;Determine the name
    
    o_path=file_dirname(out_satfix_name)
    if(get_info_stats) then begin
      ;
      n_base=strlen(out_satfix_name)
      n_dot=strpos(out_satfix_name,'.',/reverse_search)
      ;
      if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
        clog_file=strtrim(out_satfix_name,2)+'_wire_trace.log'
      endif else begin
        clog_file=strmid(strtrim(out_satfix_name,2),0,n_dot)+'_wire_trace.log'
      endelse
      print,'log file name=',clog_file
      ;
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
        print,'error opening the casing wire trace log file'
        print,'Logfile name='+strtrim(clog_file,2)
        print,'DWEL_baseline_sat_fix_wire terminating'
        print,' '
        err=7
        goto,cleanup
      endif
      ;
      printf,ctfile,strtrim('DWEL wire reflected trace Log File',2)
      printf,ctfile,strtrim('Run made at: '+systime(),2)
      printf,ctfile,strtrim('Input Cube File: '+strtrim(DWELCubeFile,2),2)
      flush,ctfile
      ;
      pos_sav=where(reform(save_w[1,*]) gt 0.0,npos_sav)
      if (npos_sav gt 0) then begin
        wire_mean=mean(reform(save_w[2,pos_sav]),/double)
        wire_stdev=stddev(reform(save_w[2,pos_sav]),/double)
        wire_cv=100.0*wire_stdev/wire_mean
        pos_mean=mean(reform(save_w[3,pos_sav]),/double)
        pos_stdev=stddev(reform(save_w[3,pos_sav]),/double)
        pos_cv=100.0*pos_stdev/pos_mean
        off_mean=mean(reform(save_w[4,pos_sav]),/double)
        off_stdev=stddev(reform(save_w[4,pos_sav]),/double)
        off_cv=100.0*off_stdev/(off_mean+0.0001)
      endif else begin
        wire_mean=0.0
        wire_stdev=0.0
        wire_cv=0.0
        pos_mean=0.0
        pos_stdev=0.0
        pos_cv=0.0
        off_mean=0.0
        off_stdev=0.0
        off_cv=0.0
      endelse
      printf,ctfile,'Number of lines used in stats='+strtrim(string(npos_sav),2)
      printf,ctfile,'Stats,Mean,CV(%)'
      outstring='Base_Stats(DN)='+strtrim(string(mean_base,format='(f14.3)'),2)+','+ $
        strtrim(string(cv_base,format='(f14.2)'),2)
      printf,ctfile,strtrim(outstring,2)
      outstring='Wire_Stats(DN)='+strtrim(string(wire_mean,format='(f14.3)'),2)+','+ $
        strtrim(string(wire_cv,format='(f14.2)'),2)
      printf,ctfile,strtrim(outstring,2)
      outstring='Pos_Stats(Bin)='+strtrim(string(pos_mean,format='(f14.3)'),2)+','+ $
        strtrim(string(pos_cv,format='(f14.2)'),2)
      printf,ctfile,strtrim(outstring,2)
      outstring='Offset_Stats(DN)='+strtrim(string(off_mean,format='(f14.3)'),2)+','+ $
        strtrim(string(off_cv,format='(f14.2)'),2)
      printf,ctfile,strtrim(outstring,2)
      outstring='wire_total(DN)='+strtrim(string(wire_power,format='(f14.2)'),2)
      printf,ctfile,strtrim(outstring,2)
      outstring='wire_fwhm(ns)='+strtrim(string(wire_fwhm,format='(f14.2)'),2)
      printf,ctfile,strtrim(outstring,2)
      printf,ctfile,'Line_Num,Wire_Num,Wire_Max_Val(DN),Wire_Max_Pos(Bin),Offset(DN),Wire_Mean_Max(DN),Wire_Power(DN),Wire_FWHM(ns),sm_Wire_Max_Val(DN)'
      flush, ctfile
      ;
      for i=0L,nl-1L do begin
        outstring=strtrim(string(save_w[0,i],format='(f14.0)'),2)+','+ $
          strtrim(string(save_w[1,i],format='(f14.0)'),2)+','+ $
          strtrim(string(save_w[2,i],format='(f14.3)'),2)+','+ $
          strtrim(string(save_w[3,i],format='(f14.3)'),2)+','+ $
          strtrim(string(save_w[4,i],format='(f14.3)'),2)+','+ $
          strtrim(string(save_w[5,i],format='(f14.3)'),2)+','+ $
          strtrim(string(save_w[6,i],format='(f14.3)'),2)+','+ $
          strtrim(string(save_w[7,i],format='(f14.3)'),2)+','+ $
          strtrim(string(sm_wire_max[i],format='(f14.3)'),2)
        ;
        printf,ctfile,outstring
      endfor
      ;
      flush, ctfile
      free_lun, ctfile,/force
      save=0b
    endif
    
    ;======================================================================
    ;now correct the casing pulse for the wire effect
    
    Twire=meantwzero/time_step
    Twire_I=round(Twire)
    x_in=float(p_time/time_step)+Twire
    x_out=float(round(p_time/time_step)+Twire_I)
    w_pulse=pulse_model*maxmeanwpulse
    help,x_in
    help,x_out
    help,w_pulse
    wp_out=interpol(w_pulse,x_in,x_out)
    x_out=round(x_out)
    pulse_corr=pulse
    pulse_corr[x_out]=pulse_corr[x_out]-wp_out
    w_model=fltarr(n_elements(pulse))
    w_model[x_out]=wp_out
    
    help,pulse_corr
    
    nu_casing_max=max(pulse_corr,nu_Tzero_I)*CasingMaxratio
    ; interpolate peak location
    istat = peak_int(wl[[nu_Tzero_I-1, nu_Tzero_I, nu_Tzero_I+1]], pulse_corr[[nu_Tzero_I-1, nu_Tzero_I, nu_Tzero_I+1]], time_int, pulse_int, offset)
    
    if (istat gt 0) then begin
      print,'peak_int had an issue for nu_Tzero from wire corrected pulse'
    endif
    
    nu_Tzero = time_int
    print, 'Initial nu_Tzero after interpolation = ', nu_Tzero, ' ns'
    
    delta= -casing2Tzero/c2 ; 0.07 meter is about the distance between the rotating mirror and the base. This is an old measurement and needs be updated.
    nu_Tzero=nu_Tzero+delta
    print,'Shifted nu_Tzero (ns)='+strtrim(string(nu_Tzero),2)
    
    nu_scale_mean=target_dn/nu_casing_max
    
    ;======================================================================
    
    ;write out the information into a spectral library to check later
    write_wire_spec_file:
    
    ;Determine library file name
    n_base=strlen(out_satfix_name)
    n_dot=strpos(out_satfix_name,'.',/reverse_search)
    if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
      wire_library=strtrim(out_satfix_name,2)+'_wire_Lib.sli'
    endif else begin
      wire_library=strmid(out_satfix_name,0,n_dot)+'_wire_Lib.sli'
    endelse
    
    ; check if we need to close and delete library file
    fids=envi_get_file_ids()
    if(fids[0] ne -1) then begin
      for i=0,n_elements(fids)-1 do begin
        envi_file_query,fids[i],fname=name
        if (strtrim(strlowcase(wire_library),2) eq $
          strtrim(strlowcase(name),2)) then begin
          envi_file_mng,id=fids[i],/remove,/delete
        endif
      endfor
    endif
    
    print,'nposw='+strtrim(string(nposw),2)
    
    out_wmat=dblarr(nposw,6)
    out_wmat=double([[wpulse],[w_model[posw]],[baseline[posw]],[wsig],[pulse[posw]],[pulse_corr[posw]]])
    out_wmat=double(reform(out_wmat,nposw,6,/overwrite))
    sense=size(out_wmat,/structure)
    help,out_wmat
    
    Lib_Nam=['wPulse','wModel','Baseline','wsig','Casing_Pulse','de_wired Casing']
    
    openw, wlun, wire_library,/get_lun
    writeu,wlun,out_wmat
    free_lun,wlun
    
    ;set up the ENVI header for the output library
    descrip='Wire Pulse Information from: '+strtrim(DWELCubeFile,2)
    envi_setup_head,fname=wire_library,ns=nposw,nl=6,nb=1,$
      xstart=0,ystart=0,file_type=4,interleave=0,$
      data_type=sense.type,spec_names=lib_nam,$
      descrip=descrip, wl=time[posw],$
      zplot_titles=['Time','Value'], $
      /write
      
    envi_open_file,wire_library,r_fid=lib_fid,/no_interactive_query,/no_realize
    envi_file_mng,id=lib_fid,/remove
    
    out_wmat=0b
  endif
  ;======================================================================
  
  DWEL_base_fix_info=strarr(21)
  DWEL_base_fix_info=[$
    'Program=DWEL_Baseline_Sat_Fix_Cmd_NSF',$
    'Descr=DWEL New Base Fix Settings with casing power',$
    'Processing Date Time='+strtrim(systime(),2),$
    'Pulse='+strtrim(pulse_model_name,2),$
    'Wire_Flag='+strtrim(string(wire,format='(i14)'),2),$
    'Comment=Tzero is the time at which the peak of the output pulse occurs',$
    'Tzero='+strtrim(string(Tzero,format='(f14.3)'),2),$
    'srate='+strtrim(string(srate,format='(f14.2)'),2),$
    'out_of_pulse='+strtrim(string(out_of_pulse,format='(i14)'),2),$
    'mean_offset='+strtrim(string(mean_line_offset,format='(f14.3)'),2),$
    'Target_dn='+strtrim(string(target_dn,format='(f14.2)'),2),$
    'scale_mean='+strtrim(string(scale_mean,format='(f14.3)'),2),$
    'Noise_RMS='+strtrim(string(mean_base_sig,format='(f14.3)'),2),$
    'Casing_Type='+strtrim(Casing_Type,2),$
    'Low(deg)='+strtrim(string(low,format='(f14.2)'),2),$
    'High(deg)='+strtrim(string(high,format='(f14.2)'),2),$
    'delta(ns)='+strtrim(string(delta,format='(f14.4)'),2),$
    'casing_MeanMax='+strtrim(string(CasingMeanWfMax,format='(f14.3)'),2),$
    'casing_MaxMean='+strtrim(string(CasingMaxWfMean,format='(f14.3)'),2),$
    'casing_Maxratio='+strtrim(string(CasingMaxratio,format='(f14.4)'),2),$
    'casing_sig='+strtrim(string(CasingMeanSig,format='(f14.3)'),2),$
    'Casing_CV(%)='+strtrim(string(CasingMeanCV,format='(f14.2)'),2),$
    'casing_fwhm(nsec)='+strtrim(string(casing_fwhm,format='(f14.4)'),2),$
    'casing_fwhm(m)='+strtrim(string(casing_fwhm*c2,format='(f14.4)'),2),$
    'model_fwhm(nsec)='+strtrim(string(model_fwhm,format='(f14.4)'),2),$
    'model_fwhm(m)='+strtrim(string(model_fwhm*c2,format='(f14.4)'),2) $
    ]
  if keyword_set(wire) then begin
    DWEL_base_fix_info=[DWEL_base_fix_info,$
      'comment='+strtrim('Wire info for correction of casing data follows:',2),$
      'wire_MeanMax='+strtrim(string(meanmaxwpulse,format='(f14.3)'),2),$
      'wire_MaxMean='+strtrim(string(maxmeanwpulse,format='(f14.3)'),2),$
      'wire_Tzero='+strtrim(string(meantwzero,format='(f14.3)'),2),$
      'corr_MaxMean='+strtrim(string(nu_casing_max,format='(f14.3)'),2),$
      'nu_scale_mean='+strtrim(string(nu_scale_mean,format='(f14.3)'),2),$
      'corr_Tzero='+strtrim(string(nu_Tzero,format='(f14.3)'),2) $
      ]
    scale_mean=nu_scale_mean
    tzero=nu_tzero
  endif
  DWEL_base_fix_info=strtrim(DWEL_base_fix_info,2)
  
  ;====================================================================
  ;now apply the correction to the image
  
  nb_out=nb
  nl_out=nl
  ns_out=ns
  wl_range=c2*time
  
  ;Open the output file for BIL tiling
  text_err=0
  
  ;====================================================================
  openw, osatfile, out_satfix_name,/get_lun,error=text_err
  if (text_err ne 0) then begin
    print, strtrim('Halting DWEL_baseline_sat_fix', 2)
    print, strtrim(['Error opening output file '+strtrim(out_satfix_name,2)], 2)
    err=8
    goto, cleanup
  endif
  ;====================================================================
  
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
  runin_mask[0:shotend-1]=0.0
  runin_mask[shotend:shottz]=merge
  runin_mask[shottz+1:nb_out-1]=1.0
  ;  help,runin_mask
  ;  mean_image=fltarr(ns_out,nl_out)
  ;  max_image=fltarr(ns_out,nl_out)
  ;====================================================================
  sat_mean_image=fltarr(ns_out,nl_out)
  sat_max_image=fltarr(ns_out,nl_out)
  satfixeddata = fltarr(ns_out,nb_out)
  ;====================================================================
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
  
  pos_stime=where(time ge -10.0 and time le 20.0,npos_stime)
  print,'number of savetimes='+strtrim(string(npos_stime),2)
  
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
    temp=float(data)
    ;==================================================================
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
    temp=temp-transpose(cmreplicate(sm_line_offset[i],nb))##one_ns
    pos_z=where(mask_all[*,i] eq 0,count_z)
    if (count_z gt 0L) then begin
      temp[pos_z,*]=0.0
    endif
    
    ;this has been done to make satfix correct as data input - ie constant baseline removed
    ;======================================================================
    satfixeddata=temp
    ;saturation fix, replace the saturated waveform part with pulse model and corresponding peak intensity
    if (count_sat gt 0) then begin
      satflag=1
      totsat=long(totsat)+long(count_sat)
      totbad=0L
      for si=0,count_sat-1 do begin
        unsattemp=reform(temp[sat_pos[si], *])
        if keyword_set(wire) then savewire=unsattemp[pos_stime]
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
          if keyword_set(wire) then unsattemp[pos_stime]=savewire
          satfixeddata[sat_pos[si], *] = unsattemp
        endelse
        unsattemp=0b
      endfor
      alltotbad=long(alltotbad)+long(totbad)
    endif
    temp=satfixeddata
    satfixeddata=0b
    
    ;======================================================================
    ;; here only use scale_mean to scale the data overall. The line-by-line
    ;; scaling is done in next fixbase stage because it is easier to deal with
    ;; wire if wire is present. 
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
    
    ;======================================================================
    ;round if integer else set back in data
    if (type lt 4 or type gt 9) then begin
      temp=fix(round(temp), type=2)
    endif
    
    ;write out the resulting tile
    writeu,osatfile,temp
    sat_mean_image[*,i]=fln*total(temp[*,pos_pos],2)
    sat_max_image[*,i]=float(max(temp[*,pos_pos],DIMENSION=2))
    ;==================================================================
    
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
  ;======================================================================
  ;get names
  out_base=file_basename(out_satfix_name)
  ;get output_ancillary file name
  ;Set up the ancillary file
  dot = strpos(out_satfix_name,'.',/reverse_search)
  if ((dot lt 0) or ((strlen(out_satfix_name)-dot-1) ne 3)) then dot = strlen(out_satfix_name)
  ancfile = strmid(out_satfix_name, 0, dot)+'_ancillary.img'
  last=strpos(ancfile,path_sep(),/reverse_search)
  anc_base=file_basename(ancfile)
  
  ;======================================================================
  ; Write saturation fixing parameters to header
  DWEL_sat_info=strarr(10)
  DWEL_sat_info = [$
    'Program=DWEL_Baseline_Sat_Fix_Cmd_NSF',$
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
  ;======================================================================
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
  ;; envi_assign_header_value, fid=anc_fid, $
  ;;   keyword='DWEL_Adaptation', $
  ;;   value=DWEL_Adaptation

  envi_assign_header_value, fid=anc_fid, keyword='DWEL_base_fix_info', $
    value=DWEL_base_fix_info
  envi_assign_header_value, fid=anc_fid, keyword='DWEL_sat_info', $
    value=DWEL_sat_info
  envi_write_file_header, anc_fid
  
  envi_file_mng,id=anc_fid,/remove
  ;====================================================================
  ;
  sat_mask=0b
  ;====================================================================
  
  ; Do the final cleanup
  cleanup:
  free_lun,lun,/force
  free_lun,inlun,/force
  free_lun,ofile,/force
  free_lun,tfile,/force
  free_lun,ctfile,/force
  free_lun,osatfile,/force
  if (err ne 0) then print,'Returning from DWEL_Baseline_Sat_Fix_Cmd_NSF with error'
  heap_gc,/verbose
  ;

  ;; write processing time summary
  print, '**************************************************'
  print, 'Processing program = dwel_baseline_sat_fix_cmd_nsf'
  print, 'Input DWEL cube file size = ' + $
    strtrim(string(double(procfilesize)/(1024.0*1024.0*1024.0)), 2) + ' G'
  print, 'Processing time = ' + strtrim(string((systime(1) - starttime)), 2) + ' ' + $
    'seconds'
  print, '**************************************************'

  return
end
