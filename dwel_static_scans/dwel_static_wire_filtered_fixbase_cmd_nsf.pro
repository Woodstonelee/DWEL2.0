;
;+
; NAME:
;DWEL_STATIC_WIRE_FILTERED_FIXBASE_CMD_NSF
;
; PURPOSE:
;This is to process a fake data cube made from stationary scan. 
;Baseline re-fix (and ancillary file update) for pre-filtered DWEL data.
;
; CATEGORY:
;DWEL waveform processing.
;
; CALLING SEQUENCE:
;dwel_static_wire_filtered_fixbase_cmd_nsf, FilteredFile, Inancfile, OutUpdatedFile, get_info_stats, zen_tweak, err
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
pro dwel_static_wire_filtered_fixbase_cmd_nsf, FilteredFile, Inancfile, OutUpdatedFile, $
  get_info_stats, zen_tweak, err, target_range

  ;; FilteredFile: the file name of the DWEL cube file that had been base and sat fixed and then filtered
  ;
  compile_opt idl2
;  envi, /restore_base_save_files
;  envi_batch_init, /no_status_window
  ;

  print, '************************************************************'
  print, 'The minimum range of stationary scans this processing, '
  print, 'dwel_static_wire_filtered_fixbase_cmd_nsf can accept = 5.0 m'
  print, 'If your stationary scan is from a closer target, '
  print, 'the processing may give unexpected results'
  print, '************************************************************'

  resolve_routine, 'DWEL_SET_THETA_PHI_NSF', /compile_full_file, /either
  resolve_routine, 'DWEL_ITPULSE_MODEL_DUAL_NSF', /compile_full_file, /either
  resolve_routine, 'DWEL_GET_HEADERS', /compile_full_file, /either
  resolve_routine, 'DWEL_PUT_HEADERS', /compile_full_file, /either
  resolve_routine, 'DT2NB', /compile_full_file, /either
  resolve_routine, 'CMREPLICATE', /compile_full_file, /either

  print, 'Target range by user: ', target_range

  ;; get the size of input file to be processed. It will be used in later
  ;; summary of processing time. 
  procfilesize = file_info(FilteredFile)
  procfilesize = procfilesize.size
  ;; get the time now as the start of processing
  starttime = systime(1)

  outrangemin = -5.0
  outrangemax = 93.5

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
    out_of_pulse=200 ; in unit of bins
  endif else begin
    out_of_pulse=400 ; in unit of bins
  endelse
  target_dn=512.0
  err=0
  err_flag=0b
  
  ;; a file to record interpolated range and intensity of wire signal (Tzero) of
  ;;all waveforms.
  wfile = 60
  wirelog_file=''
  
  ;; ===========================================================================
  ;; some internal parameters
  ;; distance from casing (edge of casing) to the true Tzero position
  ;; default is for NSF DWEL
  ;; casing2Tzero = 0.055 ; unit: meters
  ;; the FWHM of outgoing pulse, ns
  outgoing_fwhm = 5.1
  ;; the full width of outgoing pulse where intensity is below 0.01 of
  ;;maximum
  pulse_width_range = 40.0 ; in unit of ns
  print,'pulse_width_range=',pulse_width_range
  ;saturation test
  sat_test=1023L
  ;; ===========================================================================
  
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
  
  if (~status) then begin
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

  ;Read the laser manufacturer from the scan info
  ;; default laser is manlight
  match = -1
  for i=0,n_elements(DWEL_headers.DWEL_scan_info)-1 do begin
    if (strmatch(DWEL_headers.DWEL_scan_info[i],'*lasers*',/fold_case)) then match=i
  endfor
  if (match ge 0) then begin
    sf = strtrim(strcompress(strsplit(DWEL_headers.DWEL_scan_info[match],'=',/extract)),2)
    laser_man = sf[1]
  endif else begin
    laser_man = 'manlight'
  endelse  
  print,'Laser manufacturer = '+strtrim(laser_man)

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
  
  if (wavelength eq 1064) then begin
    target_dn=512.0
    ;; move leftward from out_of_pulse to find peak of wire signal
;    bins_toward_wire = 57 ; in unit of bins
  endif else begin
    target_dn=509.0
    ;; move leftward from out_of_pulse to find peak of wire signal
;    bins_toward_wire = 80 ; in unit of bins
  endelse

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

  ;update the delta/wire2Tzero from early basefix if found
  match = -1
  for i=0,n_elements(base_info)-1 do if (strmatch(base_info[i],'*delta*',/fold_case)) then match=i
  if (match ge 0) then begin
    text=strtrim(base_info[match],2)
    print,'text=',text
    k=strpos(text,'=')
    delta=float(strtrim(strmid(text,k+1),2))
    print,'info match for delta= ',strtrim(base_info[match],2)
    print,'Extracted delta=', delta
    wire2Tzero = delta*c2*(-1)
  endif

  ;; find out_of_pulse from early basefix
  match = -1
  for i=0,n_elements(base_info)-1 do if (strmatch(base_info[i],'*out_of_pulse*',/fold_case)) then match=i
  if (match ge 0) then begin
    text=strtrim(base_info[match],2)
    print,'text=',text
    k=strpos(text,'=')
    out_of_pulse=fix(strtrim(strmid(text,k+1),2))
  endif else begin

  endelse
  if (match ge 0) then print,'info match for out_of_pulse= ',strtrim(base_info[match],2)
  print,'out_of_pulse=', out_of_pulse

  
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
    ShotAzim:rotary_encoder, $
    ZenEnc:zenithenc $
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
    ;;  it's a fake data cube from stationary scan, use all waveforms in
    ;;  a scan line. 
    index=where((mask_all[*,i] ne 0) and (satmask ne 1L), count)
    if (count gt 0L) then begin
      pointsz=long64(i)*long64(bufrs)
      data=read_binary(inlun,data_start=pointsz,data_dims=[ns,nb],data_type=type)
      d = double((data[index,*]))
      data=0b
      ;; median waveform of this 'scan line'
      ;; due to severe noise in waveforms, median has been proved to be more
      ;; stable than simple mean
      temp = median(d, dimension=1, /even)
      ;; temp=total(d,1,/double)/double(count)
      ;; median of waveform bins after the wire signal, it can give us the noise
      ;;base level. 
      tmpmedian = median(d[*, out_of_pulse:nb-1], dimension=2, /even)
      ;; mean of medians of waveforms in this scan line, as the mean base of
      ;;this scan line. 
      temp2 = mean(tmpmedian, /double)
      ;; temp2=total(d[*,out_of_pulse:nb-1],/double)/(double(count)*double(nb-out_of_pulse))
      ;
      nt=n_elements(temp)
      store=reform(temp[posw])
      tempmax=max(store,nct)
      ; interpolate peak location
      istat = peak_int([float(nct-1), float(nct), float(nct+1)], store[[nct-1, nct, nct+1]], tzero_loc, value, offset)
      
      mpos=posw[0]+nct
      line_scale[i]=target_dn/(tempmax-temp2)
      t_zerol[i]=(float(posw[0])+tzero_loc)*time_step
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
;;      if ((i gt 50L) and (i lt nl-9L)) then begin
        n = long(n) + 1L
        sum = sum + temp
        sum2 = sum2 + total(d^2, 1,/double)/double(count)
;;      endif
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
  mean_base = median(pulse[out_of_pulse:nb-1],/double, /even)
  mean_base_sig = median(sig[out_of_pulse:nb-1],/double, /even)
  ;; mean_base=total(pulse[out_of_pulse:nb-1],/double)/double(nb-out_of_pulse)
  ;; mean_base_sig=total(sig[out_of_pulse:nb-1],/double)/double(nb-out_of_pulse)
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
    ;; for stationary scans, do NOT scale so that we can trace the DN for
    ;; calibration estimate
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
    printf,ctfile,'Line_Num,Casing_Num,Casing_Max_Val,Casing_Max_Pos,Offset,Line_Scale,Casing_Power,Casing_FWHM,Casing_Mean_Tzero'
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
  CasingMeanWfMax = max(reform(pulse[posw]), nct)
  Tzero_I=posw[0]+nct
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
  baseline[postarget] = 0.0
  
  pulse = pulse - baseline
  
  ;  CasingMeanWfMax = max(pulse, Tzero_I)
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
  
  ;; interpolate peak location
  istat = peak_int(time[[Tzero_I-1, Tzero_I, Tzero_I+1]], pulse[[Tzero_I-1, Tzero_I, Tzero_I+1]], time_int, pulse_int, offset)
  Tzero = float(time_int)
  CasingMeanWfMaxInt = pulse_int
  print, 'Initial Tzero after interpolation = ', Tzero, ' ns'
  
  posscal=where((t_zerol eq 0.0) or (line_scale eq 0.0),nscal)
  if (nscal gt 0) then begin
    t_zerol[posscal]=Tzero
    line_scale[posscal]=scale_mean
  endif
  
  delta= -wire2Tzero/c2 ; 0.065 meter is about the distance between the rotating mirror and the base.
  print,'delta=',delta, ' ns'
  print,''
  ;
  Tzero=Tzero+delta
  print,'Shifted Tzero='+strtrim(string(Tzero),2)+' ns'
  ;now correct the line based estimates as well
  t_zerol=t_zerol+delta
  
  ;; av_tzero=mean(t_zerol[50:nl-11])
  ;; print,'Average Tzero of each scan line = '+strtrim(string(av_tzero),2)+' ns'
  
  ;; t_zerol[0:49]=av_tzero
  ;; t_zerol[nl-10:nl-1]=av_tzero
  
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
  pulse_model_name = 'NSF_DWEL_ItPulse_Model'
  ;; if the input data is from  Oz DWEL, keopsys lasers, 
  if strcmp(laser_man, 'keopsys', /fold_case) then begin
    DWEL_itpulse_model_dual_oz, wavelength, i_val, t_val, r_val, p_range, p_time, pulse_model
    pulse_model_name = 'Oz_DWEL_ItPulse_Model'
  endif 
  print, 'Pulse model used in filtered_fixbase: '+pulse_model_name
  
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
  DWEL_filtered_fix_info=strtrim(DWEL_filtered_fix_info,2)
  
  ;=====================================================================
  ;now apply the correction to the image
  
  nb_out=nb
  nl_out=nl
  ns_out=ns
  wl_out=fltarr(nb_out)
  wl_range=c2*time

  ;; check if the available range values are within the given output range
  ;; limits between outrangemin and outrangemax. Because the digitizer start
  ;; time changes, it is possible that the digitizer starts too early or too
  ;; late that the recorded available ranges cannot cover the given output
  ;; range. If so, shrink the output range by 2 meter, more than 3 times of the
  ;; standard deviation of wire signal position obtained from stationary scans
  ;; with wire. 
  wl_range_max = max(wl_range, min=wl_range_min)
  if wl_range_max le outrangemax then begin
    outrangemax = wl_range_max - 2.0
  endif 
  if wl_range_min ge outrangemin then begin
    outrangemin = wl_range_min + 2.0
  endif 
  print, 'Output range limit: [', outrangemin, ', ', outrangemax, ']'
  
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

  runin_mask=fltarr(nb_out)+1.0
  ;; if a target range is provided, it is interpreted as the maximum possible
  ;; range of the targe in the static scan. 
  ;; suppress all waveform bins after this range to zeros to avoid complications
  ;;in later processing caused by noise. 
  if n_elements(target_range) ne 0 or arg_present(target_range) then begin
    target_time = Tzero + double(target_range)/c2
    target_end = fix((target_time + pulse_width_range)/time_step)
    target_halfend = fix((target_time + pulse_width_range/2.0)/time_step)
    runin_mask[target_end+1:nb_out-1] = 0.0
    numpos=indgen(target_end-target_halfend+1) + target_halfend
    ord=float(numpos-target_halfend)/float(target_end-target_halfend)
    merge = ((1.0-ord)^2)*(1.0+2.0*ord)
    runin_mask[target_halfend:target_end] = merge
    target_start = fix((target_time - pulse_width_range)/time_step)
    wire_end = fix((Tzero + pulse_width_range)/time_step)
    if wire_end+1 le target_start-1 then begin
      runin_mask[wire_end+1:target_start-1] = 0.0
    endif 
  endif 

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
  ;; Select a distance range for output waveforms. 
  ;; The range should be between the minimum and maximum of wl_out or even
  ;; smaller if an empirical range is supplied, e.g. [-0.6, 95.0]. 
  ;; Also notice that in the following we correct Tzero scan line by scan line
  ;; (or in future possibly shot by shot). This Tzero shift might cause the
  ;; possible maximum range in a scan line to be smaller than our select output
  ;; range here if we don't have pre-caution here. 
  ;; According to Glenn Howe in Sept. 2014, the Tzero variance from shot to shot
  ;; is around 3ns ~ 4ns. If we set the output range to be 1m within the range
  ;; of wl_out here and an empirical range, the following Tzero correction
  ;; is expected to be safe to give waveforms covering the whole output
  ;; range. Otherwise in the following processing, nb_loc (number of bins in a
  ;; waveform within output range) might not be equal to nb_resamp (number of
  ;; bins in an expected output waveform)
  
  ;; ;; first set up an empirical output range. 
  ;; outrangemin = -6.0
  ;; outrangemax = 95.0
  ;; ;; check this output range against possible range of wl_out
  ;; if outrangemin lt min(wl_out) then begin
  ;;   outrangemin = min(wl_out)
  ;; endif 
  ;; if outrangemax gt max(wl_out) then begin
  ;;   outrangemax = max(wl_out)
  ;; endif 
  ;; ;; take the Tzero variance into account
  ;; outrangemin = outrangemin + 1.0
  ;; outrangemax = outrangemax - 1.0
  ;; ;; after we got a range from a test scan, fix the output range here but do
  ;; ;; some checks and give warnings (or no warning, the following processing will
  ;; ;; throw errors and stop)
  ;; if (outrangemin - 1.0) gt -5.0 then begin
  ;;   print, 'WARNING: possible minimum range is larger than the given one!'
  ;;   print, 'Possible minimum range = ' + strtrim(string(outrangemin - 1.0), 2)
  ;;   print, 'Given minimum range = -5.0'
  ;; endif 
  ;; if (outrangemax + 1.0) lt 93.5 then begin
  ;;   print, 'WARNING: possible maximum range is smaller than the given one!'
  ;;   print, 'Possible maximum range = ' + strtrim(string(outrangemax + 1.0), 2)
  ;;   print, 'Given maximum range = 93.5'
  ;; endif 
  print, 'Output range (meter) = [' + strtrim(string(outrangemin), 2) + ', ' + $
    strtrim(string(outrangemax), 2) + ']'
  posk=where((wl_out ge outrangemin) and (wl_out le outrangemax),nb_resamp)
  print, 'Number of bins in output range = ' + strtrim(string(nb_resamp), 2)
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

  ;; see if we need to write out wire signal for log information
  if get_info_stats then begin
    n_base=strlen(OutUpdatedFile)
    n_dot=strpos(OutUpdatedFile,'.',/reverse_search)
    ;
    if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
      wirelog_file=strtrim(OutUpdatedFile,2)+'_filtfix_wire.log'
    endif else begin
      wirelog_file=strmid(strtrim(OutUpdatedFile,2),0,n_dot)+'_filtfix_wire.log'
    endelse
    print,'wire signal log file name=',wirelog_file
    ;  clog_file=o_path+path_sep()+clog_file
    ;see if the log file exists & remove if it does!
    if(file_test(wirelog_file)) then begin
      fids=envi_get_file_ids()
      if(fids[0] eq -1) then begin
      
      endif else begin
        for i=0,n_elements(fids)-1 do begin
          envi_file_query,fids[i],fname=tname
          if (strlowcase(strtrim(strlowcase(wirelog_file),2)) eq $
            strlowcase(strtrim(strlowcase(tname),2))) then begin
            envi_file_mng,id=fids[i],/remove
          endif
        endfor
        
      endelse
    endif
    ;Open wire signal log file
    text_err=0
    openw, wfile, wirelog_file,/get_lun,error=text_err
    if (text_err ne 0) then begin
      print,' '
      print,'error opening the wire signal log file'
      print,'wire signal log file name='+strtrim(wirelog_file,2)
      print,'dwel_static_wire_filtered_fixbase_cmd terminating'
      print,' '
      err_flag=1b
      err=44
      goto,cleanup
    endif
    ;
    printf,wfile,strtrim('DWEL stationary scan wire signal log File',2)
    printf, wfile, strtrim('tzero,intensity,time[k],int[k],sample,line,band,shot_num', 2)
  endif 
  
  ;do the processing over the lines for BIL structure
  shot_num = 0L
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
    temp=scale_mean*temp*(transpose(runin_mask)##one_ns)
    ;

    ;; get Tzero for each individual waveform from wire signal
    wire_tzero = fltarr(ns_out)
    lambda = fltarr(ns_out)
    ;; resampled waveform data
    retemp = fltarr(ns_out, nb_resamp)
    for j = 0, ns_out-1 do begin
      if mask_all[j,i] eq 0 then begin
        continue
      endif 
      shot_num = shot_num + 1L

      tmpwf = temp[j, *]
      tmpmax = max(tmpwf[posw], tmpmaxI)
      tmpind = posw[0] + [tmpmaxI-1, tmpmaxI, tmpmaxI+1]
      w_maxI = posw[0] + tmpmaxI
      istat = peak_int(tmpind*time_step, tmpwf[tmpind], time_int, pulse_int, $
        offset)
      wire_tzero[j] = time_int + delta

      ;; it is possible that a waveform misses wire signal and wire_tzero[j] is
      ;; not right. check wire_tzero[j] before further processing. If it is not
      ;; good, then put this shot as invalid. 
;; or (abs(pulse_int-CasingMeanWfMax) lt 3*mean_base_sig*scale_mean) or (abs(wire_tzero[j] - tzero) gt 15*time_step) 
      if (finite(time_int, /nan)) then begin
        mask_all[j, i] = 0
        wire_tzero[j] = 0
        print, 'Wire signal is bad at sample, ', j+1, ', line, ', i+1
        continue
      endif 

      ;remove wire signal from waveform
      ;get a pulse from the pulse model scaled by interpolated pulse peak
      wire_pulse = interpol(pulse_int*pulse_model, p_time, $
        p_time+tmpind[1]*time_step-time_int)
      ;subtract this wire pulse from waveform
      tmpwf[tmpind[1]-i_val[2]:tmpind[1]+i_val[n_elements(i_val)-1]-i_val[2]] = $
        tmpwf[tmpind[1]-i_val[2]:tmpind[1]+i_val[n_elements(i_val)-1]-i_val[2]] $
        - wire_pulse

      ;resample to new standard ranges
      wl_loc=fltarr(nb_out)
      wl_loc=wl_range+cmreplicate(c2*(Tzero-wire_tzero[j]),nb_out)
      posk=where(wl_loc lt 0.0,nposk)
      wlk=wl_loc[posk[nposk-1]]
      lambda[j]=-wlk/(time_step*c2)
      wl_new=wl_loc-wlk
      posk=0b
      posk=where((wl_new ge outrangemin) and (wl_new le outrangemax),nb_loc)
      if (nb_loc ne nb_resamp) then begin
        print,'Resampled shot inconsistent, nb_loc='+strtrim(string(nb_loc),2)
        print,'Expected Value='+strtrim(string(nb_resamp),2)
        print,'Scan Line Number='+strtrim(string(i+1),2)
        print, 'Shot numbeer='+strtrim(string(j+1), 2)
        print,'Tzero,wire_tzero=',Tzero,wire_tzero[j]
        print, 'Tzero intensiyt=', strtrim(string(tmpmax), 2)
        ;; now use mean tzero
        wl_loc = 0b
        wl_loc=fltarr(nb_out)
        wl_loc=wl_range
        posk=where(wl_loc lt 0.0,nposk)
        wlk=wl_loc[posk[nposk-1]]
        lambda[j]=-wlk/(time_step*c2)
        wl_new=wl_loc-wlk
        posk=0b
        posk=where((wl_new ge outrangemin) and (wl_new le outrangemax),nb_loc)
        if (nb_loc ne nb_resamp) then begin
          mask_all[j,i] = 0
          continue
        endif else begin
          wire_tzero[j] = Tzero
          pulse_int = CasingMeanWfMaxInt
          w_maxI = Tzero_I
          tmpmax = CasingMeanWfMax
        endelse 
        ;; err_flag=1b
        ;; err=33
        ;; goto,cleanup
      endif
      ; printf, wfile, strtrim('tzero,intensity,time[k],int[k],sample,line,band', 2)
      buf = string(wire_tzero[j], pulse_int, $
        (w_maxI)*time_step, tmpmax, j+1, i+1, w_maxI, $
        shot_num, $
        format='(4f14.3,4i10)') 
      buf=strtrim(strcompress(buf),2)
      while (((ii = strpos(buf, ' '))) ne -1) do $
          strput, buf, ',', ii
      printf,wfile,buf

      tmpwf = lambda * shift(tmpwf, -1) + (1.0 - lambda) * tmpwf
      retemp[j, *] = tmpwf[posk]      
    endfor 

    ;round if integer
    if (type lt 4 or type gt 9) then begin
      retemp=fix(round(retemp), type=2)
    endif
    ;write out the resulting tile
    writeu,ofile,retemp
    ;get some stats images
    mean_image[*,i]=fln*total(retemp[*,pos_pos],2)
    max_image[*,i]=float(max(retemp[*,pos_pos],DIMENSION=2))
    ;=================================================
    data=0
    temp = 0
    temp=0b
    posk=0b
    wl_loc=0b
    wire_tzero =0
    lambda = 0
    retemp = 0
  endfor
  
  flush, wfile
  free_lun, wfile, /force

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
  free_lun, wfile, /force
  if (err_flag) then print,'dwel_filtered_fixbase returned with error'
  heap_gc,/verbose
  ;
  ;; write processing time summary
  print, '**************************************************'
  print, 'Processing program = dwel_filtered_fixbase_cmd_nsf'
  print, 'Input cube file size = ' + $
    strtrim(string(double(procfilesize)/(1024.0*1024.0*1024.0)), 2) + ' G'
  print, 'Processing time = ' + strtrim(string((systime(1) - starttime)), $
    2) + ' ' + $
    'seconds'
  print, '**************************************************'

  print, '************************************************************'
  print, 'The minimum range of stationary scans this processing, '
  print, 'dwel_static_wire_filtered_fixbase_cmd_nsf can accept = 5.0 m'
  print, 'If your stationary scan is from a closer target, '
  print, 'the processing may give unexpected results'
  print, '************************************************************'

  return
end
