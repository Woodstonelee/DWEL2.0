;; baseline re-fix (and ancillary file update) for pre-filtered DWEL data
;; Zhan Li, zhanli86@bu.edu
;; Original with satfix Created in 2013 by Zhan Li
;; Last modified: 20140603 by Zhan Li
;
; This routine assumes input is filtered basefixed and satfixed file
; Inancfile will usually be the ancillary file for the basefix and satfix image (which
; is the one filtered in normal operation)
; It re-assesses the casing information and also re-writes satfix ancillary information and headers
; after this, the pulse has a new model of the Iterated Pulse and Tzero should be the position
; of the peak of the outgoing iterated pulse!
;
pro dwel_filtered_fixbase_cmd_nsf, FilteredFile, Inancfile, OutUpdatedFile, get_info_stats, zen_tweak, err, wire=wire
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
;FilteredFile = string, the file name of the DWEL cube file that had been base
;and sat fixed and then filtered with dwel_general_filter (convolution with a
;pulse model). 
;
;Inancfile = string, ancillary file name of the DWEL cube.
;
;OutUpdatedFile = string, file name of output file.
;
;get_info_stats = a boolean to ask if the program outputs statistics of the
;casing waveforms.  
;
;zen_tweak = integer, a number to tweak the zenith encoder shift.
;
;err = integer, to return error code by this program. 
;
; OUTPUTS:
;None
;
; KEYWORDS:
;wire = a binary keyword, if set, the input scan is collected with wire and the
;program will correct Tzero shot by shot with wire signal and remove wire
;signals after all processing. 
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
;David Jupp, Sept 2014, 
;  - Created this routine. 
;Zhan Li, Oct 2014,
;  - Added documentation comments.
;Zhan Li, Dec 2014, 
;  - Add one more keyword argument "wire". If it is set, it indicates the scan
;is collected with wire and will remove wire signals from waveforms after Tzero
;correction in this processing. 
;-
;

  compile_opt idl2
  ;

  ; catch IDL error caused by this routine
  catch, Error_status
  ; once we catch an error, dump out the scan line and shot number so
  ; that we can trace back directly to the bad waveform
  if Error_status ne 0 then begin
    catch, /cancel
    Help, /Last_Message, Output=theErrorMessage
    FOR j=0,N_Elements(theErrorMessage)-1 DO BEGIN
       Print, theErrorMessage[j]
    ENDFOR
    ; print, 'Error message: ', !error_state.msg
    print, 'Line, Y=', line_ind+1
    print, 'Sample, X=', sample_ind+1
    goto, cleanup
  endif 
  line_ind = -1
  sample_ind = -1

  resolve_routine, 'DWEL_SET_THETA_PHI_NSF', /compile_full_file, /either
  resolve_routine, 'DWEL_ITPULSE_MODEL_DUAL_NSF', /compile_full_file, /either
  resolve_routine, 'DWEL_GET_HEADERS', /compile_full_file, /either
  resolve_routine, 'DWEL_PUT_HEADERS', /compile_full_file, /either
  resolve_routine, 'DT2NB', /compile_full_file, /either
  resolve_routine, 'CMREPLICATE', /compile_full_file, /either
  resolve_routine, 'CLEAN_LINE', /compile_full_file, /either

  ;; get the size of input file to be processed. It will be used in later
  ;; summary of processing time. 
  procfilesize = file_info(FilteredFile)
  procfilesize = procfilesize.size
  ;; get the time now as the start of processing
  starttime = systime(1)

  lun=99
  inlun=105
  ofile=101
  tfile=98
  ctfile=35
  fname=''
  o_name=''
  ctfile=30
  before_casing=0
  check=1b
  err=0
  err_flag=0b
  wfile = 108
    
  ;; preset some values of settings, will be updated according to header
  ;; information if necessary. 
  outrangemin = -6.0
  outrangemax = 95.0

  if (check) then begin
    out_of_pulse=200
  endif else begin
    out_of_pulse=500
  endelse
  target_dn=512.0
  ; distance from casing (edge of casing) to the true Tzero position
  ;; default is the distance from scan mirror center to center of lambertian
  ;; panel at nadir point. Default here is for NSF. Will update later if
  ;; needed. 
  casing2Tzero = 0.055 ; unit: meters
  ; the FWHM of outgoing pulse, ns
  outgoing_fwhm = 5.1
  ; the full width of outgoing pulse where intensity is below 0.01 of maximum
  pulse_width_range = 40.0
  ;saturation test value in dn
  sat_test=1023L
  ;set wire limits
  w_zen_min=25.0
  w_zen_max=65.0
  ;set Filter Factor (reduction of max due to iterated filter)
  ; default is for NSF, will update later according to laser manufacturer
  FiltFactor = 0.97946 
  ;
  filt5=set_gfilt(5)
  filt10=set_gfilt(10)
  ;clean up any fids which are no longer where they were!
  ;ENVI issue that is annoying and leads to confusion
  clean_envi_file_fids
  
  print,'entering filtered basefix program'
  
  ;set speed of light metres per nsec /2
  c=0.299792458
  c2=c/2.0
  
  ; Open DWEL cube file
  envi_open_file, FilteredFile, r_fid=infile_fid,/no_realize,/no_interactive_query
  
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
  n_base=strlen(fname)
  n_dot=strpos(fname,'.',/reverse_search)
  ;  if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
  ;    anc_name=strtrim(fname,2)+'_ancillary'
  ;  endif else begin
  ;    Inancfile=strmid(fname,0,n_dot)+'_ancillary.img'
  ;  endelse
  inancfile=strtrim(inancfile)
  print,'inancfile='+inancfile
  if(~file_test(Inancfile)) then begin
    message_text=[ $
      'Ancillary file does not exist!',$
      'Expected Name='+strtrim(Inancfile,2)$
      ]
    print, message_text
    envi_file_mng,id=infile_fid,/remove
    err_flag=1b
    err=3
    GOTO, cleanup
  endif
  
  envi_open_file, Inancfile, r_fid=ancillaryfile_fid, $
    /no_realize
  ;check if operation cancelled
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

  ; update FiltFactor
  if strcmp(laser_man, 'keopsys', /fold_case) then begin
    casing2Tzero = 0.065 ; unit: meter
    FiltFactor=0.92558
  endif 
  
  Old_Tzero=0.0
  ;find the Tzero from early basefix
  match = -1
  for i=0,n_elements(base_info)-1 do if (strmatch(base_info[i],'*Tzero=*')) then match=i
  if (match ge 0) then begin
    text=strtrim(base_info[match],2)
    print,'text=',text
    k=strpos(text,'=')
    Old_Tzero=float(strtrim(strmid(text,k+1),2))
  endif else begin
    Old_Tzero=0.0
  endelse
  if (match ge 0) then print,'info match for old Tzero= ',strtrim(base_info[match],2)
  print,'Old Tzero=',Old_Tzero

  ;update the delta/casing2Tzero from early basefix if found
  match = -1
  for i=0,n_elements(base_info)-1 do if (strmatch(base_info[i],'*delta*',/fold_case)) then match=i
  if (match ge 0) then begin
    text=strtrim(base_info[match],2)
    print,'text=',text
    k=strpos(text,'=')
    delta=float(strtrim(strmid(text,k+1),2))
    print,'info match for delta= ',strtrim(base_info[match],2)
    print,'Extracted delta=', delta
    casing2Tzero = delta*c2*(-1)
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

  ;; find the casing type for laser power variation monitoring from early
  ;; basefix 
  ;; default casing type is dewired or no-wire returns from lambertian panel. 
  casing_type = 'LAM'
  match = -1
  for i=0,n_elements(base_info)-1 do if (strmatch(base_info[i],'*Casing_Type*',/fold_case)) then match=i
  if (match ge 0) then begin
    text=strtrim(base_info[match],2)
    print,'text=',text
    k=strpos(text,'=')
    casing_type=strtrim(strmid(text,k+1),2)
  endif else begin

  endelse
  if (match ge 0) then print,'info match for casing type= ',strtrim(base_info[match],2)
  print,'casing type='+casing_type

  ;; find wire flag from early basefix
  match = -1
  for i=0,n_elements(base_info)-1 do if (strmatch(base_info[i],'*Wire_Flag*',/fold_case)) then match=i
  if (match ge 0) then begin
    text=strtrim(base_info[match],2)
    print,'text=',text
    k=strpos(text,'=')
    old_wire_flag=fix(strtrim(strmid(text,k+1),2))
  endif else begin

  endelse
  if (match ge 0) then print,'info match for old wire flag= ',strtrim(base_info[match],2)
  print,'old wire flag=', old_wire_flag
  if old_wire_flag xor wire then begin
    print, 'Wire flag from previous processing is different from input wire setting here in dwel_filtered_fixbase'
    print, 'Processing still uses your input wire setting = ', wire
    print, 'Results may not be correct. Double check wire settings!'
  endif 
    
  ;set up default wl as time
  print,'initial wl[0],d_wl='+strtrim(string(wl[0]),2)+','+strtrim(string(wl[1]-wl[0]),2)
  wl=wl/c2+Old_Tzero
  wl=findgen(nb)*time_step
  print,'changed wl[0],d_wl='+strtrim(string(wl[0]),2)+','+strtrim(string(wl[1]-wl[0]),2)
  
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

  zenenc = -1
  ;; get some information for header file
  DWEL_Adaptation = ENVI_GET_HEADER_VALUE(ancillaryfile_fid, 'DWEL_Adaptation', undefined=undef)
  ;; get laser wavelength
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

    ;; find scan encoder of zenith point from header information. If not found,
    ;; then the default values in dwel_set_phi_theta will be used later.  
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
      ;; still use 512, need to manually provide calibration parameters in
      ;; dwel_get_point_cloud       
      target_dn = 512.0 
    endif 
  endif else begin
    target_dn = 509.0
    if strcmp(Casing_Type, 'CAP', /fold_case) then begin
      target_dn = 75.0
    endif
    if strcmp(Casing_Type, 'CASE', /fold_case) then begin
      ;; still use 509, need to manually provide calibration parameters in
      ;; dwel_get_point_cloud 
      target_dn = 509.0 
    endif 
  endelse
  
  if keyword_set(wire) then begin
  
    def_wire_Tzero=0.0
    ;find the mean wire Tzero from early basefix
    match = -1
    for i=0,n_elements(base_info)-1 do if (strmatch(base_info[i],'*wire_Tzero*')) then match=i
    if (match ge 0) then begin
      text=strtrim(base_info[match],2)
      k=strpos(text,'=')
      def_wire_Tzero=float(strtrim(strmid(text,k+1),2))
    endif else begin
      def_wire_Tzero=0.0
    endelse
    if (match ge 0) then print,'info match for wire_Tzero= ',strtrim(base_info[match],2)
    print,'def_wire_Tzero=',def_wire_Tzero
    
    if (def_wire_tzero eq 0.0) then print,'warning!!! no def wire tzero found!!'
    
    def_scale_mean=1.0
    ;find the Tzero from early basefix
    match = -1
    for i=0,n_elements(base_info)-1 do if (strmatch(base_info[i],'*nu_scale_mean*')) then match=i
    if (match ge 0) then begin
      text=strtrim(base_info[match],2)
      k=strpos(text,'=')
      def_scale_mean=float(strtrim(strmid(text,k+1),2))
    endif else begin
      def_scale_mean=1.0
    endelse
    if (match ge 0) then print,'info match for scale_mean= ',strtrim(base_info[match],2)
    print,'def_scale_mean=',def_scale_mean
    
    def_wire_Max=0.0
    ;find the Tzero from early basefix
    match = -1
    for i=0,n_elements(base_info)-1 do if (strmatch(base_info[i],'*wire_MaxMean*')) then match=i
    if (match ge 0) then begin
      text=strtrim(base_info[match],2)
      k=strpos(text,'=')
      def_wire_Max=float(strtrim(strmid(text,k+1),2))
    endif else begin
      def_wire_Max=0.0
    endelse
    if (match ge 0) then print,'info match for wire_Max= ',strtrim(base_info[match],2)
    print,'def_wire_Max=',def_wire_Max
    
    if (def_wire_Max eq 0.0) then print,'warning!!! no def wire Max found!!'
    
    def_wire_Max=def_wire_Max*def_scale_mean*FiltFactor
    print,'adjusted def_wire_Max=',def_wire_Max

    ;find the wire2Tzero record provided in last basefix
    ;; NOT IN USE NOW B/C TZERO IS FROM LINE-BY-LINE CASING RETURNS AFTER WIRE
    ;; SIGNAL REMOVAL RATHER THAN SHOT-BY-SHOT WIRE RETURNS.
    match = -1
    for i=0,n_elements(base_info)-1 do if (strmatch(base_info[i],'*wire2Tzero*')) then match=i
    if (match ge 0) then begin
      text=strtrim(base_info[match],2)
      k=strpos(text,'=')
      wire2Tzero=float(strtrim(strmid(text,k+1),2))
    endif else begin
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
    endelse
    if (match ge 0) then print,'info match for wire2Tzero= ',strtrim(base_info[match],2)
    print,'wire2Tzero=',wire2Tzero    
    wiredelta = wire2Tzero/c2
  endif

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
  ;
  status = dwel_set_theta_phi_nsf(p_stat,zen_tweak)
  ;put the results into the local arrays
  zeniths=(*p_stat).ShotZen
  azimuths=(*p_stat).ShotAzim
  ptr_free, p_stat
  
  ;--------------------------------------------
  ;now get the output file name
  
  output:
  
  n_base=strlen(fname)
  n_dot=strpos(fname,'.',/reverse_search)
  
  ;============================================
  
  ;check if output file is open in envi
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
  ;open input file
  openr,inlun,FilteredFile,/get_lun,error=err
  if (err gt 0) then begin
    err=7
    goto,cleanup
  endif
  
  ;=================================================================================================
  
  ;set up the pulse model
  ;; default pulse model is from NSF DWEL, manlight lasers.
  DWEL_itpulse_model_dual_nsf, wavelength, i_val, t_val, r_val, p_range, p_time, pulse_model
  pulse_model_name = 'NSF_DWEL_ItPulse_Model'
  ;; if the input data is from  Oz DWEL, keopsys lasers, 
  if strcmp(laser_man, 'keopsys', /fold_case) then begin
    DWEL_itpulse_model_dual_oz, wavelength, i_val, t_val, r_val, p_range, p_time, pulse_model
    pulse_model_name = 'Oz_DWEL_ItPulse_Model'
  endif 
  
  model_fwhm=total(pulse_model)*time_step
  
  ;======================================================================
  
  ;; get mean pulse and baseline from the given casing area designated
  ;;by the zenith angles.
  n = long(0)
  sum = dblarr(nb)
  sum2 = dblarr(nb)
  temp=dblarr(nb)
  line_scale=fltarr(nl)+1.0
  t_zerol=fltarr(nl)
  casing_max=fltarr(nl)
  casing_num=lonarr(nl)
  casing_offset=fltarr(nl)
  if (get_info_stats) then begin
    save=fltarr(9,nl)
    if keyword_set(wire) then save_w=fltarr(8,nl)
  endif
  if keyword_set(wire) then begin
    nw=long(0)
    w_start=old_tzero-0.4-pulse_width_range/4.0
    w_end=old_tzero-0.4+pulse_width_range/4.0
    posw=where((wl ge w_start) and (wl le w_end),nposw)
    sumw = dblarr(nposw)
    sumw2 = dblarr(nposw)
    save_tpos=fltarr(nl)
    save_tmax=fltarr(nl)
    save_valid=fltarr(nl)
  endif
  count=long(0)
  old_mean_base=0.0
  ;set up the pointer for read_binary
  bufrs=long64(nbytes)*long64(nb)*long64(ns)
  pointsz=long64(0)
  ;
  for i=0L,nl-1L do BEGIN
    ;set the satmask
    satmask=long(reform(anc_data[*,i,0]))
    ;read in the data for line i
    pointsz=long64(i)*long64(bufrs)
    data=read_binary(inlun,data_start=pointsz,data_dims=[ns,nb],data_type=type)
    
    ;if there is a wire get wire stats
    if keyword_set(wire) then begin
    
      index=where((mask_all[*,i] gt 0) and ((zeniths[*,i] ge w_zen_min) and (zeniths[*,i] LE w_zen_max) $
        and (satmask eq 0b)), count)
      if (count le 50) then begin
        count=0L
        index=0b
        goto,wire_nodata
      endif
      if (count gt 0) then begin
        ;
        d = double(data[index,*])
        ;      d=reform(d[*,posw])
        nw = long(nw) + 1L
        temp=total(d,1,/double)/double(count)
        ;
        ;      nt=n_elements(temp)
        tempmax=max(temp[posw],nct)
        ;new section
        ;print,'first tempmax,nct=',tempmax,nct
        ; new code to get peak nearest centre
        r=(temp[posw]-old_mean_base)/(tempmax-old_mean_base)
        test_pos=[replicate(0.0,nposw),r,replicate(0.0,nposw)]
        ;      filt=set_gfilt(2)
        ;      test_pos=convol(test_pos,filt)
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
        tempmax=temp[posw[nct]]
        ;now see if there are two adjacent likely max positions and select max one
        if (nump gt 1) then begin
          if (npt eq 0) then begin
            if (peaks[npt+1] eq peaks[npt]+1) then begin
              if (temp[posw[peaks[npt+1]]] gt tempmax) then npt=npt+1
            endif
          endif else if (npt eq nump-1) then begin
            if (peaks[npt-1] eq peaks[npt]-1) then begin
              if (temp[posw[peaks[npt-1]]] gt tempmax) then npt=npt-1
            endif
          endif else begin
            if (peaks[npt+1] eq peaks[npt]+1) then begin
              if (temp[posw[peaks[npt+1]]] gt tempmax) then npt=npt+1
            endif else if (peaks[npt-1] eq peaks[npt]-1) then begin
              if (temp[posw[peaks[npt-1]]] gt tempmax) then npt=npt-1
            endif
          endelse
        endif
        nct=peaks[npt]
        tempmax=temp[posw[nct]]
        ;print,'output tempmax,nct=',tempmax,nct
        r=0b
        obj=0b
        ;back from new section
        mpos=posw[nct]
        truepos=float(mpos)
        ; interpolate peak location
        istat = peak_int([float(mpos-1), float(mpos), float(mpos+1)], temp[[mpos-1, mpos, mpos+1]], truepos, value, offset)
        tempmax=(1.0d0*temp[mpos-1]+2.0d0*temp[mpos]+1.0d0*temp[mpos+1])/4.0d0
        ;
        save_tpos[i]=float(truepos)*time_step
        save_tmax[i]=float(tempmax-old_mean_base)
        save_valid[i]=1b
        if (get_info_stats) then begin
          save_w[0,i]=float(i)
          save_w[1,i]=float(count)
          save_w[2,i]=save_tmax[i]
          save_w[3,i]=save_tpos[i]/time_step
          save_w[4,i]=float(old_mean_base)
          save_w[5,i]=save_tmax[i]
          save_w[6,i]=float((total((temp[posw]-old_mean_base),/double)))
          save_w[7,i]=(wl[1]-wl[0])*save_w[6,i]/save_w[2,i]
        endif
        sumw = sumw + temp[posw]
        sumw2 = sumw2 + total(d[*,posw]^2, 1,/double)/double(count)
      endif else begin
        wire_nodata:
        save_valid[i]=0b
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
    endif
    
    ;de-allocate
    index=0b
    temp=0b
    
    ;now for casing stats
    index=where((mask_all[*,i] ne 0) and ((zeniths[*,i] ge low) and (zeniths[*,i] LE high)) $
      and (satmask ne 1L), count)
    if (count ge 50L) then begin
      ;
      d = double(reform(data[index,*]))
      temp=total(d,1,/double)/double(count)
      temp2=total(d[*,out_of_pulse:nb-1],/double)/(double(count)*double(nb-out_of_pulse))
      nt=n_elements(temp)
      
      ;now remove the wire pulse if wire present
      if keyword_set(wire) then begin
        if (save_valid[i]) then begin
          set_wire_Tzero=save_tpos[i]
          set_wire_Max=save_tmax[i]
        endif else begin
          set_wire_Tzero=def_wire_Tzero
          set_wire_Max=def_wire_Max
        endelse
        Twire=set_wire_Tzero/time_step
        Twire_I=round(Twire)
        x_in=float(p_time/time_step)+Twire
        x_out=float(round(p_time/time_step)+Twire_I)
        w_pulse=pulse_model*set_wire_Max
        wp_out=interpol(w_pulse,x_in,x_out)
        x_out=round(x_out)
        temp[x_out]=temp[x_out]-wp_out
      endif
      
      store=reform(temp[before_casing:nt-1])
      tempmax=max(store,nct)
      ; interpolate peak location
      istat = peak_int([float(nct-1), float(nct), float(nct+1)], store[[nct-1, nct, nct+1]], tzero_loc, value, offset)
      
      mpos=before_casing+nct
      casing_num[i]=count
      ;; the casing_max here is after wire being subtracted if wire is
      ;; presented. i.e. dewired or no-wire casing_max
      casing_max[i]=tempmax-temp2
      line_scale[i]=target_dn/(tempmax-temp2)
      t_zerol[i]=(float(before_casing)+tzero_loc)*time_step
      casing_offset[i]=temp2
      if (get_info_stats) then begin
        save[0,i]=float(i)
        save[1,i]=float(count)
        save[2,i]=tempmax-temp2
        save[3,i]=float(mpos)
        save[4,i]=float(temp2)
        save[5,i]=float(line_scale[i])
        save[6,i]=float((total((temp-temp2),/double)))
        save[7,i]=(wl[1]-wl[0])*save[6,i]/save[2,i]
        save[8,i]=t_zerol[i]
      endif
      if ((i gt 50L) and (i lt nl-10L)) then begin
        n = long(n) + 1L
        sum = sum + temp
        sum2 = sum2 + total(d^2, 1,/double)/double(count)
      endif
    endif else begin
      count=0L
      index=0b
      casing_num[i]=0
      casing_max[i]=0.0
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
    ;
    index=0b
    data=0b
    d=0b
    temp=0b
    temp2=0b
    store=0b
  endfor
  
  ;make some space
  d=0b
  data=0b
  temp=0b
  temp2=0b
  index=0b
  
  ;now process the casing stats
  d=double(n)
  pulse=float(sum/d)
  sig=float(sqrt(abs(sum2-sum^2/d)/d))
  ;help,sig
  ;make more space
  sum=0b
  sum2=0b
  
  ;mu and sig are length the number of bands
  mean_base=total(pulse[out_of_pulse:nb-1],/double)/double(nb-out_of_pulse)
  mean_base_sig=total(sig[out_of_pulse:nb-1],/double)/double(nb-out_of_pulse)
  if (abs(mean_base) lt 1.0e-3) then cv_base=0.0 else cv_base=100.0*mean_base_sig/mean_base
  
  pulse=pulse-mean_base
  ;casing_power=sqrt(total(pulse^2,/double))
  casing_power=total(pulse[0:out_of_pulse],/double)
  if (abs(casing_power) lt 1.0e-6) then casing_fwhm=0.0 else $
    casing_fwhm=(wl[1]-wl[0])*casing_power/max(pulse[0:out_of_pulse])
  ;
    
  ;Now form the smoothed data and see what it is like!
  ;
  valid_all=bytarr(nl)+1b
  
  if keyword_set(wire) then begin
    ;first the wire if present
    wire_valid=valid_all
    sm_tpos=fltarr(nl)
    sm_tmax=fltarr(nl)
    posv=where((save_tpos le 0.0) or (save_tmax le 0.0) or ~save_valid,nposv)
    if (nposv gt 0) then wire_valid[posv]=0b
    clean_line, save_tpos,wire_valid,sm_tpos,filt5,err
    clean_line, save_tmax,wire_valid,sm_tmax,filt5,err
    help,sm_tpos
    help,sm_tmax
    posv=0b
  endif
  ;
  ;now the casing
  casing_valid=valid_all
  sm_casing_max=fltarr(nl)
  sm_t_zerol=fltarr(nl)
  temp=fltarr(nl)
  posv=where((t_zerol le 0.0) or (casing_max le 0.0) or (casing_num le 0),nposv)
  if (nposv gt 0) then casing_valid[posv]=0b
  clean_line, t_zerol,casing_valid,temp,filt5,err
  clean_line, temp,valid_all,sm_t_zerol,filt10,err
  ;
  sm_casing_offset=fltarr(nl)
  clean_line,casing_offset,casing_valid,sm_casing_offset,filt5,err
  ;
  temp=fltarr(nl)
  clean_line, casing_max,casing_valid,temp,filt5,err
  clean_line, temp,valid_all,sm_casing_max,filt10,err
  help,sm_t_zerol
  help,sm_casing_max
  posv=0b
  temp=0b
  ;
  ;; the casing_max here is after wire being subtracted if wire is
  ;; presented. i.e. dewired or no-wire casing_max 
  sm_line_scale=replicate(target_dn,nl)/sm_casing_max
  
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
      print,'DWEL_Filtered_FixBase_Cmd terminating'
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
      pos_cv=100.0*pos_stdev/pos_mean
      off_mean=mean(reform(save[4,pos_sav]),/double)
      off_stdev=stddev(reform(save[4,pos_sav]),/double)
      off_cv=100.0*off_stdev/(off_mean+0.0001)
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
    outstring='Mean Casing_Stdev(DN)='+strtrim(string(casing_power,format='(f14.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    outstring='Casing_fwhm(ns)='+strtrim(string(casing_fwhm,format='(f14.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    printf,ctfile,'Line_Num,Casing_Num,Casing_Max_Val(DN),Casing_Max_Pos(Bin),Offset(DN),Line_Scale,Casing_Power(DN),Casing_FWHM(ns),Tzero_L(ns),sm_Cmax(DN),sm_Cpos(ns),sm_Line_Scale'
    flush, ctfile
    ;
    for i=0L,nl-1L do begin
      outstring=strtrim(string(save[0,i],format='(f14.0)'),2)+','+ $
        strtrim(string(save[1,i],format='(f14.0)'),2)+','+ $
        strtrim(string(save[2,i],format='(f14.3)'),2)+','+ $
        strtrim(string(save[3,i],format='(f14.0)'),2)+','+ $
        strtrim(string(save[4,i],format='(f14.3)'),2)+','+ $
        strtrim(string(save[5,i],format='(f14.3)'),2)+','+ $
        strtrim(string(save[6,i],format='(f14.3)'),2)+','+ $
        strtrim(string(save[7,i],format='(f14.3)'),2)+','+ $
        strtrim(string(save[8,i],format='(f14.3)'),2)+','+ $
        strtrim(string(sm_casing_max[i],format='(f14.3)'),2)+','+ $
        strtrim(string(sm_t_zerol[i],format='(f14.3)'),2)+','+ $
        strtrim(string(sm_line_scale[i],format='(f14.3)'),2)
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
  ;  CasingMeanWfMax = max(reform(pulse[before_casing:n_elements(pulse)-1]), nct)
  CasingMaxWfMean=max(reform(pulse[before_casing:n_elements(pulse)-1]), nct)
  posnl=where(casing_num gt 0L,nposnl)
  CasingMeanWfMax=mean(casing_max[posnl],/double)
  casingMaxratio=CasingMeanWfMax/CasingMaxWfMean
  Tzero_I=before_casing+nct
  print, 'Initial Tzero before baseline fix = ', time[Tzero_I], ' ns'
  tlow=time[Tzero_I] - 1.5*pulse_width_range
  thigh=time[Tzero_I] + 1.5*pulse_width_range
  print,'tlow,thigh=',tlow,thigh
  baseline = dblarr(nb)
  tmpind = where(time LT tlow, count)
  IF count GT 0 THEN BEGIN
    baseline[tmpind] = pulse[tmpind]
  ENDIF
  tmpind = where(time GT thigh, count)
  IF count GT 0 THEN BEGIN
    baseline[tmpind] = pulse[tmpind]
  ENDIF
  tmpind=where((time ge tlow) and (time le thigh),count)
  if (count gt 0) then begin
    baseline[tmpind]=0.0
  endif
  
  pulse = pulse - baseline
  
  ;  CasingMeanWfMax = max(pulse, Tzero_I)
  CasingMeanSig=sig[Tzero_I]
  CasingMeanCV=100.0*CasingMeanSig/CasingMeanWfMax
  
  casing_power=total(reform(pulse[0:out_of_pulse]),/double)
  if (abs(casing_power) lt 1.0e-6) then casing_fwhm=0.0 else $
    casing_fwhm=(wl[1]-wl[0])*casing_power/max(pulse[0:out_of_pulse])
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
  sm_t_zerol=sm_t_zerol+delta
  
  av_tzero=mean(sm_t_zerol[51:nl-11])
  print,'Average Tzero='+strtrim(string(av_tzero),2)+' ns'
  t_zero=av_tzero
  ;  t_zerol[0:50]=av_tzero
  ;  t_zerol[nl-10:nl-1]=av_tzero
  
  time=time-av_Tzero
  
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
  
  
  ;now for the wire
  ;
  if keyword_set(wire) then begin
  
    if (nw le 0) then begin
      print,'No valid casing values at all!!'
      err=67
      goto,cleanup
    endif
    dw=double(nw)
    wpulse=float(sumw/dw)-mean_base
    wsig=float(sqrt(abs(sumw2-sumw^2/dw)/dw))
    
    sumw=0b
    sumw2=0b
    
    maxmeanwpulse=max(wpulse,kmax)
    
    posmaxw=where(save_valid gt 0b,nposmaxw)
    meanmaxwpulse=mean(save_tmax[posmaxw],/double)
    meantwzero=mean(save_tpos[posmaxw],/double)
    ;
    wire_power=total(wpulse,/double)*time_step
    if (abs(wire_power) lt 1.0e-6) then wire_fwhm=0.0 else $
      wire_fwhm=(wl[1]-wl[0])*wire_power/maxmeanwpulse
      
    print,'wire_fwhm='+strtrim(wire_fwhm,2)
    
    ;if get_info_stats set then write out the saved data to stats file
    ;Determine the name
    
    o_path=file_dirname(OutUpdatedFile)
    if(get_info_stats) then begin
      ;
      n_base=strlen(OutUpdatedFile)
      n_dot=strpos(OutUpdatedFile,'.',/reverse_search)
      ;
      if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
        clog_file=strtrim(OutUpdatedFile,2)+'_wire_trace.log'
      endif else begin
        clog_file=strmid(strtrim(OutUpdatedFile,2),0,n_dot)+'_wire_trace.log'
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
      printf,ctfile,strtrim('Input Cube File: '+strtrim(FilteredFile,2),2)
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
      printf,ctfile,'Line_Num,Wire_Num,Wire_Max_Val(DN),Wire_Max_Pos(Bin),Offset(DN),Wire_Mean_Max(DN),Wire_Power(DN),Wire_FWHM(ns),sm_Wmax(DN),sm_Wpos(ns)'
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
          strtrim(string(sm_tmax[i],format='(f14.3)'),2)+','+ $
          strtrim(string(sm_tpos[i],format='(f14.3)'),2)
        ;
        printf,ctfile,outstring
      endfor
      ;
      flush, ctfile
      free_lun, ctfile,/force
      save=0b
    endif
    
    ;===========================================================
    ;now get the model wire pulse to plot
    ;when correction working by line remove pulse_corr
    Twire=meantwzero/time_step
    Twire_I=round(Twire)
    x_in=float(p_time/time_step)+Twire
    x_out=float(round(p_time/time_step)+Twire_I)
    w_pulse=pulse_model*maxmeanwpulse
    wp_out=interpol(w_pulse,x_in,x_out)
    x_out=round(x_out)
    ;pulse_corr=pulse
    ;pulse_corr[x_out]=pulse_corr[x_out]-wp_out
    w_model=fltarr(n_elements(pulse))
    w_model[x_out]=wp_out
    
    ;nu values will not be needed after correction for wire by line
    ;nu_casing_max=max(pulse_corr,nu_Tzero_I)*CasingMaxratio
    ; interpolate peak location
    ;istat = peak_int(wl[[nu_Tzero_I-1, nu_Tzero_I, nu_Tzero_I+1]], pulse_corr[[nu_Tzero_I-1, nu_Tzero_I, nu_Tzero_I+1]], time_int, pulse_int, offset)
    ;if (istat gt 0) then begin
    ;  print,'peak_int had an issue for nu_Tzero from wire corrected pulse'
    ;endif
    ;nu_Tzero = time_int
    ;print, 'Initial nu_Tzero after interpolation = ', nu_Tzero, ' ns'
    ;delta= -casing2Tzero/c2 ; 0.07 meter is about the distance between the rotating mirror and the base. This is an old measurement and needs be updated.
    ;nu_Tzero=nu_Tzero+delta
    ;print,'Shifted nu_Tzero (ns)='+strtrim(string(nu_Tzero),2)
    ;nu_scale_mean=target_dn/nu_casing_max
    ;==================================================================
    ;
    ;write out the information into a spectral library to check later
    write_wire_spec_file:
    
    ;Determine library file name
    n_base=strlen(OutUpdatedFile)
    n_dot=strpos(OutUpdatedFile,'.',/reverse_search)
    if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
      wire_library=strtrim(OutUpdatedFile,2)+'_wire_Lib.sli'
    endif else begin
      wire_library=strmid(OutUpdatedFile,0,n_dot)+'_wire_Lib.sli'
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
    
    out_wmat=dblarr(nposw,5)
    out_wmat=double([[wpulse],[w_model[posw]],[baseline[posw]],[wsig],[pulse[posw]]])
    out_wmat=double(reform(out_wmat,nposw,5,/overwrite))
    sense=size(out_wmat,/structure)
    help,out_wmat
    
    Lib_Nam=['wPulse','wModel','Baseline','wsig','Casing_Pulse']
    
    openw, wlun, wire_library,/get_lun
    writeu,wlun,out_wmat
    free_lun,wlun
    
    ;set up the ENVI header for the output library
    descrip='Wire Pulse Information from: '+strtrim(FilteredFile,2)
    envi_setup_head,fname=wire_library,ns=nposw,nl=5,nb=1,$
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
  
  DWEL_filtered_fix_info=strarr(21)
  DWEL_filtered_fix_info=[$
    'Program=dwel_filtered_fixbase_cmd_nsf',$
    'Descr=DWEL New Filtered Base Fix Settings',$
    'Processing Date Time='+strtrim(systime(),2),$
    'Pulse='+strtrim(pulse_model_name, 2),$
    'Wire_Flag='+strtrim(string(wire,format='(i14)'),2),$
    'Comment=Tzero is the time at which the peak of the output iterated pulse occurs',$
    'Tzero='+strtrim(string(Tzero,format='(f14.3)'),2),$
    'srate='+strtrim(string(srate,format='(f14.2)'),2),$
    'out_of_pulse='+strtrim(string(out_of_pulse,format='(i14)'),2),$
    'Target_dn='+strtrim(string(target_dn,format='(f14.2)'),2),$
    'scale_mean='+strtrim(string(scale_mean,format='(f14.3)'),2),$
    'Noise_RMS='+strtrim(string(mean_base_sig,format='(f14.3)'),2),$
    'Casing_Type='+strtrim(Casing_Type,2),$
    'Low(deg)='+strtrim(string(low,format='(f14.2)'),2),$
    'High(deg)='+strtrim(string(high,format='(f14.2)'),2),$
    'delta(ns)='+strtrim(string(delta,format='(f14.4)'),2),$
    'casing_max='+strtrim(string(CasingMeanWfMax,format='(f14.3)'),2),$
    'casing_sig='+strtrim(string(CasingMeanSig,format='(f14.3)'),2),$
    'Casing_CV(%)='+strtrim(string(CasingMeanCV,format='(f14.2)'),2),$
    'casing_fwhm(nsec)='+strtrim(string(casing_fwhm,format='(f14.4)'),2),$
    'casing_fwhm(m)='+strtrim(string(casing_fwhm*c2,format='(f14.4)'),2),$
    'model_fwhm(nsec)='+strtrim(string(model_fwhm,format='(f14.4)'),2),$
    'model_fwhm(m)='+strtrim(string(model_fwhm*c2,format='(f14.4)'),2) $
    ]
    
  if keyword_set(wire) then begin
    DWEL_filtered_fix_info=[DWEL_filtered_fix_info,$
      'comment='+strtrim('Wire info for correction of casing data follows (wire_Tzero NOT delta corrected):',2),$
      'wiredelta(ns)='+strtrim(string(wiredelta,format='(f14.3)'),2),$
      'wire_MeanMax='+strtrim(string(meanmaxwpulse,format='(f14.3)'),2),$
      'wire_MaxMean='+strtrim(string(maxmeanwpulse,format='(f14.3)'),2),$
      'wire_Tzero(ns)='+strtrim(string(meantwzero,format='(f14.3)'),2) $
      ;      'corr_MaxMean='+strtrim(string(nu_casing_max,format='(f14.3)'),2),$
      ;      'nu_scale_mean='+strtrim(string(nu_scale_mean,format='(f14.3)'),2),$
      ;      'corr_Tzero='+strtrim(string(nu_Tzero,format='(f14.3)'),2) $
      ]
  ;    scale_mean=nu_scale_mean
  ;    tzero=nu_tzero
  endif
  
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
    print, strtrim('Halting dwel_filtered_fixbase_cmd', 2)
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
  ;
  ;==========================================
  mean_image=fltarr(ns_out,nl_out)
  max_image=fltarr(ns_out,nl_out)
  temp = fltarr(ns_out,nb_out)
  ;==========================================
  ;
  wl_out=fltarr(nb_out)
  posk=where(wl_range lt 0.0,nposk)
  wlk=wl_range[posk[nposk-1]]
  lambda=-wlk/(time_step*c2)
  wl_out=wl_range-wlk
  posk=0b
  posk=where((wl_out ge outrangemin) and (wl_out le outrangemax),nb_resamp)
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

  if keyword_set(wire) then begin
    wire_max_time = fltarr(ns_out,nl_out)
    wire_max = fltarr(ns_out, nl_out)
  endif

  ;; a counter to record Tzero from wire that causes trouble, possibly due to
  ;; its location too far away from the average Tzero or a odd second pulse
  ;; close to it, or simply due to the wire signal is missing from this shot. 
  nbadwiretzero = 0
  
  for i=0, nl_out-1 do begin
    line_ind = i
    ;first get the data tile
    data=read_binary(inlun,data_start=pointsz,data_dims=[ns_out,nb_out],data_type=type)
    pointsz=long64(pointsz)+long64(bufrs)
    pos_z=where(mask_all[*,i] eq 0,count_z)
    temp=float(data)
    temp=temp-transpose(baseline+sm_casing_offset[i])##one_ns
    if (count_z gt 0L) then begin
      temp[pos_z,*]=0.0
    endif    
    ;
    ;remove wire pulse if wire present
    if keyword_set(wire) then begin
      ;; ***********************************************************************
      ;; correct Tzero line by line
      set_wire_Tzero=sm_tpos[i]
      set_wire_Max=sm_tmax[i]
      Twire=set_wire_Tzero/time_step
      Twire_I=round(Twire)
      x_in=float(p_time/time_step)+Twire
      x_out=float(round(p_time/time_step)+Twire_I)
      w_pulse=pulse_model*set_wire_Max
      wp_out=interpol(w_pulse,x_in,x_out)
      x_out=round(x_out)
      temp[*,x_out]=temp[*,x_out]-wp_out##one_ns
      if (count_z gt 0L) then begin
        temp[pos_z,*]=0.0
      endif    
      ;scale the data to standard power level
      temp=sm_line_scale[i]*temp
      ;resample to new standard ranges
      wl_loc=fltarr(nb_out)
      wl_loc=wl_range+cmreplicate(c2*(av_Tzero-sm_t_zerol[i]),nb_out)
      posk=where(wl_loc lt 0.0,nposk)
      wlk=wl_loc[posk[nposk-1]]
      lambda=-wlk/(time_step*c2)
      wl_new=wl_loc-wlk
      posk=0b
      posk=where((wl_new ge outrangemin) and (wl_new le outrangemax),nb_loc)
      if (nb_loc ne nb_resamp) then begin
        print,'Resampled shot inconsisten, nb_loc='+strtrim(string(nb_loc),2)
        print,'Expected Value='+strtrim(string(nb_resamp),2)
        print,'IDL Line Number='+strtrim(string(i),2)
        print,'Tzero,T_zerol=',av_Tzero,sm_t_zerol[i]
        err_flag=1b
        err=33
        goto,cleanup
      endif
      ;
      temp=lambda*shift(temp,0,-1)+(1.0-lambda)*temp
      temp=reform(temp[*,posk])
      ;; end of correct Tzero line by line
      ;; ***********************************************************************

      
      ;; ;; ***********************************************************************
      ;; ;; correct Tzero shot by shot.
      ;; ;; resampled waveform data in a scan line
      ;; retemp = fltarr(ns_out, nb_resamp)
      ;; w_start_bin=fix((sm_tpos[i]-0.4-pulse_width_range/4.0)/time_step)
      ;; w_end_bin=fix((sm_tpos[i]-0.4+pulse_width_range/4.0)/time_step)
      ;; for j = 0, ns_out-1 do begin
      ;;   sample_ind = j
      ;;   if mask_all[j,i] eq 0 then begin
      ;;     continue
      ;;   endif 
      ;;   tmpmax = max(temp[j, w_start_bin:w_end_bin], tmpmaxI)
      ;;   ;; in very rare cases, a waveform could miss the wire signal
      ;;   ;; pulse
      ;;   if tmpmax gt 3*mean_base_sig then begin
      ;;     tmpind = w_start_bin + [tmpmaxI-1, tmpmaxI, tmpmaxI+1]
      ;;     istat = peak_int(tmpind*time_step, temp[j, tmpind], time_int, pulse_int, $
      ;;       offset)
      ;;     wire_max_time[j, i] = time_int
      ;;     wire_max[j, i] = pulse_int
      ;;     ;remove wire signal from waveform
      ;;     ;get a pulse from the pulse model scaled by interpolated pulse peak
      ;;     set_wire_Tzero=time_int
      ;;     set_wire_Max=sm_tmax[i]
      ;;     Twire=set_wire_Tzero/time_step
      ;;     Twire_I=round(Twire)
      ;;     x_in=float(p_time/time_step)+Twire
      ;;     x_out=float(round(p_time/time_step)+Twire_I)
      ;;     w_pulse=pulse_model*set_wire_Max
      ;;     wp_out=interpol(w_pulse,x_in,x_out)
      ;;     x_out=round(x_out)
      ;;     goodind = where(x_out ge 0, tmpcount)
      ;;     temp[j, x_out[goodind]] = temp[j, x_out[goodind]] - wp_out[goodind]          
      ;;     wl_loc=fltarr(nb_out)
      ;;     wl_loc=wl_range+cmreplicate(c2*(av_Tzero-set_wire_Tzero-wiredelta),nb_out)
      ;;   endif else begin
      ;;     print, 'Warning, detected intensity is too low to be accepted as wire signal'
      ;;     print, 'Line (X)=', i+1
      ;;     print, 'Sample (Y)=', j+1
      ;;     print, 'Detected wire intensity=', tmpmax
      ;;     wire_max_time[j, i]=0
      ;;     wire_max[j, i]=0
      ;;     wl_loc=fltarr(nb_out)
      ;;     wl_loc=wl_range+cmreplicate(c2*(av_Tzero-sm_t_zerol[i]),nb_out)
      ;;     nbadwiretzero = nbadwiretzero + 1
      ;;   endelse 

      ;;   posk=where(wl_loc lt 0.0,nposk)
      ;;   wlk=wl_loc[posk[nposk-1]]
      ;;   lambda=-wlk/(time_step*c2)
      ;;   wl_new=wl_loc-wlk
      ;;   posk=0b
      ;;   posk=where((wl_new ge outrangemin) and (wl_new le outrangemax),nb_loc)
      ;;   if (nb_loc ne nb_resamp) then begin
      ;;     print,'Resampled shot inconsistent, nb_loc='+strtrim(string(nb_loc),2)
      ;;     print,'Expected Value='+strtrim(string(nb_resamp),2)
      ;;     print,'Scan Line Number='+strtrim(string(i+1),2)
      ;;     print, 'Shot numbeer='+strtrim(string(j+1), 2)
      ;;     print,'Tzero,wire_max_time=',av_Tzero,wire_max_time[j, i]
      ;;     print, 'wire_max=', strtrim(string(wire_max[j, i]), 2)
      ;;     ;; now use average casing line Tzero if last attemp uses wire Tzero
      ;;     if wire_max_time[j, i] ne 0 then begin
      ;;       wl_loc = 0b
      ;;       wl_loc=fltarr(nb_out)
      ;;       wl_loc=wl_range+cmreplicate(c2*(av_Tzero-sm_t_zerol[i]),nb_out)
      ;;       ;resample to new standard ranges
      ;;       posk=where(wl_loc lt 0.0,nposk)
      ;;       wlk=wl_loc[posk[nposk-1]]
      ;;       lambda=-wlk/(time_step*c2)
      ;;       wl_new=wl_loc-wlk
      ;;       posk=0b
      ;;       posk=where((wl_new ge outrangemin) and (wl_new le outrangemax),nb_loc)
      ;;       nbadwiretzero = nbadwiretzero + 1
      ;;     endif 
      ;;     ;; now check the nb_loc and nb_resamp again
      ;;     if (nb_loc ne nb_resamp) then begin
      ;;       print, 'Average casing line Tzero failed resampling'
      ;;       ;; mask_all[j,i] = 0
      ;;       ;;continue
      ;;       err_flag=1b
      ;;       err=33
      ;;       goto,cleanup
      ;;     endif 
      ;;   endif
      ;;   ;
      ;;   temp[j, *]=sm_line_scale[i]*temp[j, *]
      ;;   temp[j, *]=lambda*shift(temp[j, *],0,-1)+(1.0-lambda)*temp[j, *]
      ;;   retemp[j, *] = temp[j, posk]
      ;;   wl_loc = 0b
      ;; endfor 
      ;; temp = retemp
      ;; ;; end of correct Tzero shot by shot.
      ;; ;; ***********************************************************************
    endif else begin
      ;scale the data to standard power level
      ;
      temp=sm_line_scale[i]*temp
      ;
      ;resample to new standard ranges
      wl_loc=fltarr(nb_out)
      wl_loc=wl_range+cmreplicate(c2*(av_Tzero-sm_t_zerol[i]),nb_out)
      posk=where(wl_loc lt 0.0,nposk)
      wlk=wl_loc[posk[nposk-1]]
      lambda=-wlk/(time_step*c2)
      wl_new=wl_loc-wlk
      posk=0b
      posk=where((wl_new ge outrangemin) and (wl_new le outrangemax),nb_loc)
      if (nb_loc ne nb_resamp) then begin
        print,'Resampled shot inconsisten, nb_loc='+strtrim(string(nb_loc),2)
        print,'Expected Value='+strtrim(string(nb_resamp),2)
        print,'IDL Line Number='+strtrim(string(i),2)
        print,'Tzero,T_zerol=',av_Tzero,sm_t_zerol[i]
        err_flag=1b
        err=33
        goto,cleanup
      endif
      ;
      temp=lambda*shift(temp,0,-1)+(1.0-lambda)*temp
      temp=reform(temp[*,posk])
    endelse 

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

  print, 'Number of detected bad wire returns if wire is used = ', nbadwiretzero

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
  ;anc_data=0b
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
  
  ;if info_stats is set write out a final action to check the result of the normalisations
  if (get_info_stats) then begin
    ;re-open the output file
    err=0
    openr,inlun,OutUpdatedFile,/get_lun,error=err
    if (err gt 0) then begin
      err=17
      goto,cleanup
    endif
    ;; get mean pulse and baseline from the given casing area designated
    ;;by the zenith angles.
    n = long(0)
    sum = dblarr(nb_resamp)
    sum2 = dblarr(nb_resamp)
    temp=dblarr(nb_resamp)
    line_scale=fltarr(nl)+1.0
    t_zerol=fltarr(nl)
    casing_max=fltarr(nl)
    casing_num=lonarr(nl)
    save=fltarr(9,nl)
    casing_max_pos_sd = fltarr(nl)
    count=long(0)
    ;set up the pointer for read_binary
    bufrs=long64(nbytes)*long64(nb_resamp)*long64(ns)
    pointsz=long64(0)
    ;
    for i=0L,nl-1L do BEGIN
      ;read in the data for line i
      pointsz=long64(i)*long64(bufrs)
      data=read_binary(inlun,data_start=pointsz,data_dims=[ns,nb_resamp],data_type=out_type)
      ;; get mean pulse and baseline from the given casing area designated
      ;;by the zenith angles.
      ;casing stats
      index=where((mask_all[*,i] ne 0) and ((zeniths[*,i] ge low) and (zeniths[*,i] LE high)), count)
      if (count ge 50L) then begin
        ;
        d = double(reform(data[index,*]))
        ;; check the variation of casing max pos in a scan line
        all_casing_max = max(d, all_casing_max_ind, dimension=2)
        temp = array_indices(d, all_casing_max_ind)
        all_casing_max_pos = reform(temp[1, *])
        casing_max_pos_sd[i] = stddev(all_casing_max_pos)*time_step
        temp = 0b
        ;; end of checking variation of casing max pos in a scan line
        temp=total(d,1,/double)/double(count)
        temp2=total(d[*,out_of_pulse:nb_resamp-1],/double)/(double(count)*double(nb_resamp-out_of_pulse))
        nt=n_elements(temp)
        ;
        store=reform(temp[before_casing:nt-1])
        tempmax=max(store,nct)
        ; interpolate peak location
        istat = peak_int([float(nct-1), float(nct), float(nct+1)], store[[nct-1, nct, nct+1]], tzero_loc, value, offset)
        mpos=before_casing+nct
        casing_num[i]=count
        casing_max[i]=tempmax-temp2
        line_scale[i]=target_dn/(tempmax-temp2)
        t_zerol[i]=(float(before_casing)+tzero_loc)*time_step
        save[0,i]=float(i)
        save[1,i]=float(count)
        save[2,i]=tempmax-temp2
        save[3,i]=float(mpos)
        save[4,i]=float(temp2)
        save[5,i]=float(line_scale[i])
        save[6,i]=float((total((temp-temp2),/double)))
        save[7,i]=(wl[1]-wl[0])*save[6,i]/save[2,i]
        save[8,i]=t_zerol[i]
        if ((i gt 50L) and (i lt nl-9L)) then begin
          n = long(n) + 1L
          sum = sum + temp
          sum2 = sum2 + total(d^2, 1,/double)/double(count)
        endif
      endif else begin
        count=0L
        index=0b
        casing_num[i]=0
        casing_max[i]=0.0
        line_scale[i]=0.0
        t_zerol[i]=0.0
        save[0,i]=float(i)
        save[1,i]=0.0
        save[2,i]=0.0
        save[3,i]=0.0
        save[4,i]=0.0
        save[5,i]=1.0
        save[6,i]=0.0
        save[7,i]=0.0
        save[8,i]=0.0
      endelse
      ;
      index=0b
      data=0b
      d=0b
      temp=0b
      store=0b
    endfor
    
    ;make some space
    d=0b
    data=0b
    temp=0b
    temp2=0b
    index=0b
    free_lun,inlun,/force
    ;now process the casing stats
    d=double(n)
    pulse=float(sum/d)
    sig=float(sqrt(abs(sum2-sum^2/d)/d))
    ;help,sig
    ;make more space
    sum=0b
    sum2=0b
    
    ;mu and sig are length the number of bands
    mean_base=total(pulse[out_of_pulse:nb_resamp-1],/double)/double(nb_resamp-out_of_pulse)
    mean_base_sig=total(sig[out_of_pulse:nb_resamp-1],/double)/double(nb_resamp-out_of_pulse)
    if (abs(mean_base) lt 1.0e-3) then cv_base=0.0 else cv_base=100.0*mean_base_sig/mean_base
    
    pulse=pulse-mean_base
    ;casing_power=sqrt(total(pulse^2,/double))
    casing_power=total(pulse[0:out_of_pulse],/double)
    if (abs(casing_power) lt 1.0e-6) then casing_fwhm=0.0 else $
      casing_fwhm=(wl[1]-wl[0])*casing_power/max(pulse[0:out_of_pulse])
    ;
    ;Now form the smoothed casing data and see what it is like!
    casing_valid=valid_all
    sm_casing_max=fltarr(nl)
    sm_t_zerol=fltarr(nl)
    temp=fltarr(nl)
    posv=where((t_zerol le 0.0) or (casing_max le 0.0) or (casing_num le 0),nposv)
    if (nposv gt 0) then casing_valid[posv]=0b
    clean_line, t_zerol,casing_valid,temp,filt5,err
    clean_line, temp,valid_all,sm_t_zerol,filt10,err
    temp=fltarr(nl)
    clean_line, casing_max,casing_valid,temp,filt5,err
    clean_line, temp,valid_all,sm_casing_max,filt10,err
    help,sm_t_zerol
    help,sm_casing_max
    posv=0b
    temp=0b
    ;
    sm_line_scale=replicate(target_dn,nl)/sm_casing_max
    
    ;=======================================================
    ;write out the saved data
    
    o_path=file_dirname(OutUpdatedFile)
    ;strtrim(FilteredFile,2)
    n_base=strlen(OutUpdatedFile)
    n_dot=strpos(OutUpdatedFile,'.',/reverse_search)
    ;
    if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
      clog_file=strtrim(OutUpdatedFile,2)+'_filtfix_final_check_trace.log'
    endif else begin
      clog_file=strmid(strtrim(OutUpdatedFile,2),0,n_dot)+'_filtfix_final_check_trace.log'
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
      print,'DWEL_Filtered_FixBase_Cmd terminating'
      print,' '
      err_flag=1b
      err=9
      goto,cleanup
    endif
    ;
    printf,ctfile,strtrim('DWEL calibration Output Summary casing trace Log File',2)
    printf,ctfile,strtrim('Run made at: '+systime(),2)
    printf,ctfile,strtrim('Cube File: '+strtrim(OutUpdatedFile,2),2)
    flush,ctfile
    ;
    pos_sav=where(reform(save[1,*]) gt 0.0,npos_sav)
    if (npos_sav gt 0) then begin
      casing_mean=mean(reform(save[2,pos_sav]),/double)
      casing_stdev=stddev(reform(save[2,pos_sav]),/double)
      casing_cv=100.0*casing_stdev/casing_mean
      pos_mean=mean(reform(save[3,pos_sav]),/double)
      pos_stdev=stddev(reform(save[3,pos_sav]),/double)
      pos_cv=100.0*pos_stdev/pos_mean
      off_mean=mean(reform(save[4,pos_sav]),/double)
      off_stdev=stddev(reform(save[4,pos_sav]),/double)
      off_cv=100.0*off_stdev/(off_mean+0.0001)
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
    outstring='Scale_Stats(DN)='+strtrim(string(old_scale_mean,format='(f14.3)'),2)+','+ $
      strtrim(string(old_scale_cv,format='(f14.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    outstring='sm_Scale_Stats(DN)='+strtrim(string(scale_mean,format='(f14.3)'),2)+','+ $
      strtrim(string(scale_cv,format='(f14.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    outstring='Mean Casing_Stdev(DN)='+strtrim(string(casing_power,format='(f14.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    outstring='Casing_fwhm(ns)='+strtrim(string(casing_fwhm,format='(f14.2)'),2)
    printf,ctfile,strtrim(outstring,2)
    printf,ctfile,'Line_Num,Casing_Num,Casing_Max_Val(DN),Casing_Max_Pos(Bin),Offset(DN),Line_Scale,Casing_Power(DN),Casing_FWHM(ns),Tzero_L(ns),sm_CMax(DN),sm_Cpos(ns),Casing_Max_Pos_SD(ns),sm_Line_scale'
    flush, ctfile
    ;
    for i=0L,nl-1L do begin
      outstring=strtrim(string(save[0,i],format='(f14.0)'),2)+','+ $
        strtrim(string(save[1,i],format='(f14.0)'),2)+','+ $
        strtrim(string(save[2,i],format='(f14.3)'),2)+','+ $
        strtrim(string(save[3,i],format='(f14.0)'),2)+','+ $
        strtrim(string(save[4,i],format='(f14.3)'),2)+','+ $
        strtrim(string(save[5,i],format='(f14.3)'),2)+','+ $
        strtrim(string(save[6,i],format='(f14.3)'),2)+','+ $
        strtrim(string(save[7,i],format='(f14.3)'),2)+','+ $
        strtrim(string(save[8,i],format='(f14.3)'),2)+','+ $
        strtrim(string(sm_casing_max[i],format='(f14.3)'),2)+','+ $
        strtrim(string(sm_t_zerol[i],format='(f14.3)'),2)+','+ $
        strtrim(string(casing_max_pos_sd[i],format='(f14.3)'),2)+','+ $
        strtrim(string(sm_line_scale[i],format='(f14.3)'),2)
      ;
      printf,ctfile,outstring
    endfor
    ;
    flush, ctfile
    free_lun, ctfile,/force
    save=0b
    out_mat=0b
    baseline_corr=0b
    pulse=0b
  endif
  ;===========================================================

  ;===========================================================
  ;; if a scan with wire and asking for stats information output, write wire
  ;; tzero and intensity to an image file.  
  if get_info_stats and keyword_set(wire) then begin
    dot = strpos(OutUpdatedFile,'.',/reverse_search)
    if ((dot lt 0) or ((strlen(OutUpdatedFile)-dot-1) ne 3)) then dot = strlen(OutUpdatedFile)
    wirefile = strmid(OutUpdatedFile, 0, dot)+'_wire.img'
    wire_base=file_basename(wirefile)

    ;; write wire image file.
    text_err=0
    openw, wfile, wirefile,/get_lun,error=text_err
    if (text_err ne 0) then begin
      print, strtrim('Halting dwel_filtered_fixbase', 2)
      print, strtrim(['Error opening output file '+strtrim(wire_base,2)], 2)
      err_flag=1b
      err=14
      goto, cleanup
    endif

    for j=0,3 do begin
      writeu,wfile,long(anc_data[*,*,j])
    endfor
    writeu,wfile,round(100.0*wire_max_time)
    writeu,wfile,round(wire_max)
    writeu,wfile,long(mask_all)
    writeu,wfile,round(100.0*zeniths)
    writeu,wfile,round(100.0*azimuths)

    ENVI_SETUP_HEAD, fname=wirefile, $
      ns=ns_out, nl=nl_out, nb=9, $
      interleave=0, data_type=3, $
      /write, $
      bnames=['Sat Mask','Sun Mask','Scan Encoder','Rotary Encoder', $
      'WireMax Time ns*100','WireCasing Peak','Mask','Zenith','Azimuth']

    envi_open_file,wirefile,r_fid=wire_fid,/no_interactive_query,/no_realize

    ;write out the previous header records
    status=DWEL_put_headers(wire_fid,DWEL_headers)
    ;
    ;write the new header(s) into the HDR file
    envi_assign_header_value, fid=wire_fid, keyword='DWEL_filtered_fix_info', $
      value=DWEL_filtered_fix_info
    envi_write_file_header, wire_fid

    envi_file_mng,id=wire_fid,/remove

  endif 
  ;===========================================================

  ; Final release of memory
  anc_data = 0b

  ;=====================================================================
  ; Do the final cleanup
  
  cleanup:
  free_lun, lun,/force
  free_lun, inlun,/force
  free_lun, ofile,/force
  free_lun, tfile,/force
  free_lun, ctfile,/force
  free_lun, wfile, /force
  if (err_flag) then print,'dwel_filtered_fixbase_cmd returned with error'
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

  return
end
