;
;+
; NAME:
;DWEL_SWOP_PULSE_XC_NSF
;
; PURPOSE:
;Correlate the pulse model with waveforms and then fix the Tzero-shift after
;cross correlation. In such way each return pulse will become an iterated pulse
;with symmetrical side lobes and some noise may be filtered. But the transient
;ringing noise could be enhanced. 
;
; CATEGORY:
;DWEL waveform processing.
;
; CALLING SEQUENCE:
;dwel_swop_pulse_xc_nsf, inbsfixfile, inbsfixancfile, outxcfile, zen_tweak, ierr
;
; INPUTS:
;inbsfixfile = input file name of baseline and saturation fixed waveform cube. 
;inbsfixancfile = input file name of baseline and saturation fixed ancillary
;data. 
;outxcfile = output file name of the cross-correlation results. 
;zen_tweak = a value in encoder unit to tweak the zenith point of encoder
;measurements. 
;ierr = a variable to receive error code.
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
;First, cross-correlate waveforms with a pulse model of the OTHER wavelength,
;i.e. correlate 1064 waveform data with 1548 pulse model and vice versa. 
;Second, fix Tzero-shift due to cross correlation with an asymmetric pulse
;model.  
;
; MODIFICATION HISTORY:
;David Jupp, Sept 2014 - Created a predecessor to this routine. 
;Zhan Li, Oct 2014 - Created this routine based on David Jupp's old routine.
;-
pro dwel_swop_pulse_xc_nsf, inbsfixfile, inbsfixancfile, outxcfile, zen_tweak, $
  ierr, wire=wire

  compile_opt idl2
;  envi, /restore_base_save_files
;  envi_batch_init, /no_status_window

  resolve_routine, 'DWEL_GENERAL_FILTER', /compile_full_file, /either
  resolve_routine, 'DWEL_FILTERED_FIXBASE_CMD_NSF', /compile_full_file, /either
  resolve_routine, 'DWEL_GET_HEADERS', /compile_full_file, /either
  resolve_routine, 'DWEL_PULSE_MODEL_DUAL_NSF', /compile_full_file, /either

  err_flag = 0

  ;; get the size of input file to be processed. It will be used in later
  ;; summary of processing time. 
  procfilesize = file_info(inbsfixfile)
  procfilesize = procfilesize.size
  ;; get the time now as the start of processing
  starttime = systime(1)
  
  BasefixImageFile = inbsfixfile
  AncillaryFile = inbsfixancfile
  outcube_filter = outxcfile
  
  ;Open DWEL base_sat_fixed cube file
  envi_open_file, BasefixImageFile, r_fid=infile_fid, /no_realize, $
    /no_interactive_query
  if (infile_fid eq -1) then begin
    print,strtrim('Error opening input file = '+strtrim(BasefixImageFile),2)
    err_flag=1b
    goto,out
  endif
  envi_file_query, infile_fid, ns=ns, nl=nl, nb=nb, wl=wl, $
    xstart=xstart, ystart=ystart, data_type=type, $
    interleave=ftype, fname=fname, dims=dims
    
  f_base = file_basename(fname)
  ;; set up a base structure for the DWEL headers
  DWEL_headers = { $
    f_base:f_base $
    }
  ;; find all of the DWEL headers in the hdr file
  status = dwel_get_headers(infile_fid, DWEL_headers)
  if (not status) then begin
    print,strtrim('Bad call to DWEL_get_headers!!',2)
    print,'File='+strtrim(fname,2)
    envi_file_mng,id=infile_fid,/remove
    err_flag=1b
    goto,out
  endif
  if (DWEL_headers.headers_present le 0s or not DWEL_headers.run_present) then begin
    print,strtrim('Input file is NOT a DWEL file',2)
    print,'File='+strtrim(fname,2)
    envi_file_mng,id=infile_fid,/remove
    err_flag=1b
    goto,out
  endif
  info=DWEL_headers.dwel_adaptation
  ;; now get the DWEL wavelength
  match = -1
  for i=0,n_elements(info)-1 do begin
    if (strmatch(info[i],'*Wavelength=*', /fold_case)) then match=i
  endfor
  if match ge 0 then begin
    text=strtrim(info[match],2)
    kpos=strpos(text,'=')
    wavelength=fix(strtrim(strmid(text,kpos+1,4),2))
  endif else begin
    print,strtrim('Input file does NOT have a DWEL file wavelength',2)
    print,'File='+strtrim(fname,2)
    envi_file_mng,id=infile_fid,/remove
    err_flag=1b
    goto,out
  endelse
  envi_file_mng,id=infile_fid,/remove
  wl_string=strtrim(string(wavelength),2)
  if (strpos(fname,wl_string) lt 0) then begin
    print,'Warning - Wavelength from file header does not appear in file name'
    print,'Wavelength from Header='+wl_string
    print,'File='+strtrim(fname,2)
  endif

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
  
  ;swop the filters so that data are filtered with other wavelength pulse
  sel_wl=wavelength
  if (wavelength eq 1064) then sel_wl=1548 else sel_wl=1064
  
  ;set up the pulse model
  ;; default pulse model is from NSF DWEL, manlight lasers.
  DWEL_pulse_model_dual_nsf, sel_wl, i_val, t_val, r_val, p_range, p_time, pulse, t_fwhm, r_fwhm
  pulse_model_name = 'NSF_DWEL_Pulse_Model'
  ;; if the input data is from  Oz DWEL, keopsys lasers, 
  if strcmp(laser_man, 'keopsys', /fold_case) then begin
    DWEL_pulse_model_dual_oz, sel_wl, i_val, t_val, r_val, p_range, p_time, pulse, t_fwhm, r_fwhm
    pulse_model_name = 'Oz_DWEL_Pulse_Model'
  endif 

  print,''
  print,'Number of values in filter='+strtrim(string(n_elements(pulse)),2)
  
  pulse=pulse/total(pulse)
  ierr=0
  maxval=max(pulse,mpos)
  
  print,'Mpos='+strtrim(string(mpos),2)
  ;run a general filter (not dwel specific)
  istat = dwel_general_filter(BasefixImageFile, pulse, mpos, outcube_filter, $
    ierr)
    
  if(~file_test(outcube_filter) or (ierr ne 0) or (istat ne 0)) then begin
    print,'Output file from DWEL_general_filter does NOT exist or error!'
    print,'expected file is='+strtrim(outcube_filter,2)
    if (ierr ne 0) then print,'Local error code='+strtrim(string(ierr),2)
    err_flag=1b
    goto,out
  endif
  
  print,'Output Filtered Image='+strtrim(outcube_filter,2)
  
  ;now update all the stats and ancillary file etc
  ;new routine to finalise all the calibrations, Tzero etc
  ;; set up the file name of the updated cube
  outupdatedfile = strtrim(strmid(outcube_filter,0,strpos(outcube_filter,'.', $
    /reverse_search))+'_update.img',2)
  ;; note ancillary file is from previous run before dwel_general_filter but
  ;; file is filtered by the general filter
  ierr = 0
  get_info_stats = 1
  print, 'Start re-fixing cross-correlation results and write update file ...' 

  if keyword_set(wire) then begin
    dwel_filtered_fixbase_cmd_nsf, outcube_filter, AncillaryFile, outupdatedfile, $
      get_info_stats, zen_tweak, ierr, /wire
  endif else begin
    dwel_filtered_fixbase_cmd_nsf, outcube_filter, AncillaryFile, outupdatedfile, $
      get_info_stats, zen_tweak, ierr
  endelse 
    
  if (ierr gt 0) then begin
    print,'DWEL_Filtered_FixBase_Cmd returned with error'
    print,'Local error number='+strtrim(string(ierr),2)
    err_flag=1b
    goto,out
  endif
  
  if(~file_test(outupdatedfile)) then begin
    print,'Output updated file does NOT exist or error!'
    print,'expected file is='+strtrim(outupdatedfile,2)
    err_flag=1b
    goto,out
  endif
  
  print,'Output updated Filtered Image='+strtrim(outupdatedfile,2)
  
  ;Set up expected ancillary file
  
  AncillaryFile = strtrim(strmid(outupdatedfile,0,strpos(outupdatedfile,'.', $
    /reverse_search))+'_ancillary.img',2)
  if(~file_test(AncillaryFile)) then begin
    print,'Expected Output ancillary file from DWEL_filtered_fixbase_cmd_nsf ' + $
      'does NOT exist!'
    print,'expected file is='+strtrim(AncillaryFile,2)
    goto,out
    err_flag=1b
  endif
  ;
  
  print,'Output ancillary updated Filtered Image file created: ' + $
    ''+strtrim(AncillaryFile,2)
    
  out:
  if (err_flag) then begin
    print,'An ERROR has occurred! Check the IDL console outputs'
  endif
  ierr = err_flag
  heap_gc,/verbose
  
  ;; write processing time summary
  print, '*******************************************'
  print, 'Processing program = dwel_swop_pulse_xc_nsf'
  print, 'Input cube file size = ' + $
    strtrim(string(double(procfilesize)/(1024.0*1024.0*1024.0)), 2) + ' G'
  print, 'Processing time = ' + strtrim(string((systime(1) - starttime)), $
    2) + ' ' + $
    'seconds'
  print, '*******************************************'

end
