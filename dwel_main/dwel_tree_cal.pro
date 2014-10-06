;======================================================================
;  The following set of procedures is part of the DWEL processing chain.
;  They follow the structure of the CSIRO EOC Hyperion Workshop module set.
;  The routines all have a similar structure.
;  They can be used stand-alone or as part of a Project structure with
;  an ENVI main menu button.
;
;  The general structure for a procedure called "name" is:
;
;  [Procedures called from name-doit]
;    These may be called anything
;    they are arranged so that every procedure occurs before its first
;    reference. This is good IDL practice. They all use the compile_opt idl2
;    convention.
;  pro name_doit
;    This is the main area of action for the algorithm - the DOIT!
;  [Procedures supporting the general top level procedure]
;      These procedures are not part of the algorithm - just the framework.
;      The workshop routines have a common top level look and feel.
;    pro name_banner_text
;      Top level widget Help and how to use the module
;    pro name_disclaimer
;      General disclaimer to all responsibility and blame accessed
;      from main widget
;    pro name_action_go
;      Main widget event handler
;    pro name_action_exit
;      Another event handler
;    pro name_resize
;      Yet another event handler
;    pro name_cleanup
;      Will they never end?
;    pro name
;      This is the common top level procedure that simply sets up an
;      IDL widget structure that allows the user to look at the main
;      help and if desired hit GO and start the module. It handles
;      exit nicely.
;======================================================================
;
function set_tree_cal_params, pstat
compile_opt idl2

def_thresh=(*pstat).fwhm_threshold
def_scale=(*pstat).mean_tree_DBH
def_zero_hit=(*pstat).save_rois
def_F=(*pstat).F
def_cal=(*pstat).calibrate
def_refl=(*pstat).tree_refl

set_base:

;first set up the widget base and user help
wb_fix_base=widget_auto_base( $
  title='Set Tree Calibrate Parameters',$
  group=(*pstat).top)

Info_Text=[$
'The following information is provided for the tree calibration routine',$
' ',$
'FWHM Threshold (m) [2.5] is a threshold defining FWHM for a hard hit.',$
'The default value is conservatively larger than expected for a',$
'hard hit. Shots that have larger FWHM than this level will be removed from',$
'the average in a ROI',$
' ',$
'Mean Tree DBH (m) [0.3] is an estimate for the mean DBH of the trees in',$
'the stand from which the sample trees are selected. This is useed to allow',$
'for trunk curvature in the reflectance model for the tree returns.',$
'Default is 0.3m',$
' ',$
'Calibration scale factor (F) is the zero range laser power factor obtained',$
'from the casing values in the new version of basefix. If the factor is not',$
'in the headers it is set to 1.0 and should be entered manually. It brings',$
'the mean casing power to a standard value and allows for changes in',$
'outgoing laser power or system transmittance.',$
' ',$
'Tree reflectance [0.4] is the estimated diffuse reflectance of a tree trunk',$
'assuming normal laser incidence and a flat surface. The default is based on',$
'field data but should be replaced with measured values if available.',$
' ',$
'Save ROI means [0] determines if the mean values of the return over the ROI',$
'are saved as a spectral library. This allows the mean return to be examined',$
'more carefully if there seems to be a problem with the data.',$
'',$
'Calibrate [0] is a flag to undertake the calibration factors estimation. It',$
'is usual to first run with Save ROI Means set and check the results carefully.',$
'When all is well run with No Save ROI means and Calibrate set. The calibrate',$
'results are recorded in an ASCII file. THis file normally has to be further',$
'processed in a spreadsheet since some errors are hard to filter out. The file',$
'is a comma separated values file.',$
' ' $
]
;put in the help and some next level bases to use
;for input info
w_fix_txt=widget_slabel(wb_fix_base,prompt=Info_Text,$
  /frame,xsize=50,ysize=10)

;set threshold

wb_L1=widget_base(wb_fix_base,/column,/frame)

wb_left=widget_base(wb_L1,/row,/frame)

wb_l2=widget_param(wb_left,/auto_manage,default=def_thresh,$
         Prompt='FWHM Threshold (m): ', $
         uvalue='fwhm_threshold',floor=0.0,ceil=50.0,field=4)

wb_l3=widget_param(wb_left,/auto_manage,default=def_scale,$
         Prompt='Mean tree DBH (m): ', $
         uvalue='mean_tree_DBH',floor=0.0,ceil=10.0,field=4)

wb_mid=widget_base(wb_L1,/row,/frame)

wb_l4=widget_param(wb_mid,/auto_manage,default=def_F,$
Prompt='Cal Scale Factor (F)',$
uvalue='F',floor=0.0,ceil=100.0,field=4)

wb_l4=widget_param(wb_mid,/auto_manage,default=def_refl,$
Prompt='Tree Reflectance: ',$
uvalue='tree_refl',floor=0.0,ceil=1.0,field=3)

wb_right=widget_base(wb_L1,/column,/frame)

wb_first=widget_base(wb_right,/row)

list=['Save ROI means in SL','Do not Save ROI means']

;set the save zero hits option
w_signature=widget_menu(wb_first,/auto_manage,default_ptr=def_zero_hit,$
         Prompt='Save ROIs: ',list=list,/exclusive,$
         uvalue='save_rois',rows=6)

wb_next=widget_base(wb_L1,/column,/frame)

list=['Calibrate','Do NOT Calibrate']

w_evi=widget_menu(wb_first,/auto_manage,default_ptr=def_cal,$
         Prompt='Calibrate Option: ',list=list,/exclusive,$
         uvalue='calibrate',rows=6)

;realise the widget using ENVI auto-manage
result=auto_wid_mng(wb_fix_base)

;result is a structure - see if setup cancelled
if (result.accept eq 0) then return,0b

def_thresh=result.fwhm_threshold
def_F=result.F

(*pstat).fwhm_threshold=result.fwhm_threshold
(*pstat).mean_tree_DBH=result.mean_tree_DBH
(*pstat).save_rois=result.save_rois
(*pstat).tree_refl=result.tree_refl
(*pstat).F=result.F
(*pstat).calibrate=result.calibrate

return,1b

end

function phase_trunk,width,dbh,range,r_crit
compile_opt idl2
;
; Function for the trunk phase effect due to curvature
; assume width is in mrad, dbh in metres, range in metres

p2=width/2000.0
if(p2 lt 1.0e-6) then begin
  phase_trunk=1.0
  return,phase_trunk
endif
r_crit=dbh*(1.0-sin(p2))/(2.0*sin(p2))
if (range gt r_crit) then begin
  print,'range is outside r_crit'
  range_loc=r_crit
endif else range_loc=range

th2=asin(2.0*(range_loc+dbh/2.0)*sin(p2)/dbh)
phase_trunk=sin(th2-p2)/(th2-p2)

return,phase_trunk
end

pro DWEL_tree_cal_doit, event
; Filter DWEL data based on pulse shape
compile_opt idl2

tfile=98
outfile=''
o_name=''
ofile=99L

;set up catch
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
   result=dialog_message(info_text,/error,title='Error in DWEL_tree_cal_doit')
  goto, cleanup
endif

;clean up any fids which are no longer where they were!
;ENVI issue that is annoying and leads to confusion
clean_envi_file_fids

;now clean up any rois around on the system!
roi_ids = envi_get_roi_ids()
if (roi_ids[0] ge 0)then begin
  print,roi_ids[0],' current roi ids found, deleted'
  envi_delete_rois,/all
endif

infile=''
anc_name=''

print,'starting dwel tree cal'

input:
; Select the input file - either from within ENVI or by opening the file
infile_select:
envi_select, /file_only, fid=fid, file_type=0,/no_dims, /no_spec, $
  title='Select DWEL data file for generating a tree calibration'
if (fid eq -1) then begin
    result=dialog_message('Try again ? (No/Yes) ',/question,$
    title='No file selected or operation cancelled !')
    if(result eq 'No') then goto, cleanup else goto, infile_select
endif

; Get the file dimensions etc
envi_file_query, fid, fname=infile, nb=nbands, nl=nlines, $
  ns=nsamples, bnames=bnames, wl=range, data_type=dt

num_range=n_elements(range)

n_base=strlen(infile)
n_dot=strpos(infile,'.',/reverse_search)
f_path = file_dirname(infile)

;get image name as separate string
last=strpos(infile,path_sep(),/reverse_search)
f_base=strtrim(strmid(infile,last+1,n_base-last-1),2)

if ((num_range le 1) or (num_range ne nbands)) then begin
  print,'num_range=',num_range
  print,'nbands=',nbands
  info_text=[ $
  'File: '+strtrim(f_base),$
  ' may not be an DWEL Laser Data file!!',$
  '',$
  'The number of ranges is not equal to the number of bands',$
  'or there is only one range in this file',$
  'It may be an Ancillary data file - check it!',$
  'Attempting to use these data for calibration will not work.',$
  '',$
  'Really use this file for tree calibration ? (No/Yes) ' $
  ]
  result=dialog_message(info_text,/question,$
  title='File may not be an DWEL Data file!')
  if(result eq 'No') then goto, infile_select
endif

;set up a base structure for the DWEL headers
  DWEL_headers={ $
     f_base:f_base $
     }
;find all of the DWEL headers in the hdr file as defined by FID
  status=dwel_get_headers(fid,DWEL_headers)

;respond to major problems!
  if (~status) then begin
    result=dialog_message('Bad FID for selected file - Try again ? (Yes.No)',/question,$
      title='bad FID in DWEL get_headers on input file')
    envi_file_mng,id=fid,/remove
    if (result eq 'Yes') then begin
      goto, infile_select
    endif else begin
      goto, cleanup
    endelse
  endif
  if (DWEL_headers.headers_present le 0s or DWEL_headers.run_present le 0s) then begin
    result=dialog_message('Selected file invalid DWEL file - Try again ? (Yes.No)',/question,$
      title='Input File is NOT a valid DWEL file')
    envi_file_mng,id=fid,/remove
    if (result eq 'Yes') then begin
      goto, infile_select
    endif else begin
      goto, cleanup
    endelse
  endif

;check for processing Level
level1='File has been baseline fixed'
level2='File has been saturation fixed'
if (~(DWEL_headers.base_present) or ~(DWEL_headers.sat_present)) then begin
  if (~(DWEL_headers.base_present)) then level1='File has NOT been baseline fixed'
  if (~(DWEL_headers.sat_present)) then level2='File has NOT been saturation fixed'
  info_text=[ $
  'File header indicates DWEL data may not have been processed to the',$
  'recommended Level for Tree Calibration. The recommended Level is for',$
  'DWEL data to be baseline fixed, saturation fixed and filtered.',$
  'For the File '+f_base+':',$
  level1,level2,$
  '',$
  'Really use this file for tree calibration ? (No/Yes) ' $
  ]
  result=dialog_message(info_text,/question,$
  title='File has NOT been processed to recommended level')
  if(result eq 'No') then goto, infile_select
endif

got_cm=0b

if (~(DWEL_headers.base_present)) then begin
  CM=132.77
  scale=1.0
endif else begin
; Get the casing max from the Base_fix info.
  info = DWEL_headers.DWEL_base_fix_info
  match = -1
  for i=0,n_elements(info)-1 do if (strmatch(info[i],'*casing_max*')) then match=i
  if (match ge 0) then begin
    sf = strsplit(info[match],'=',/extract)
    CM = float(sf[1])
    got_cm=1b
  endif else begin
    print,'casing_max not found in dwel header'
    CM=132.77
  endelse
  F=132.77/CM
;
print,'casing max=',CM
print,'scale_change=',F

if (~got_cm) then begin
info_text=[$
'The Casing Maximum not found in headers',$
'File may be an old version or not correct',$
'Try a new input file ? (Yes/No)']
  result=dialog_message(info_text,/question,$
    title='Casing Maximum not found in headers')
  envi_file_mng,id=fid,/remove
  if (result eq 'Yes') then begin
    goto, infile_select
  endif else begin
    goto, cleanup
  endelse
endif

got_scale=0b

;now get the data scale for the base fixing
  info = DWEL_headers.DWEL_base_fix_info
  match = -1
  for i=0,n_elements(info)-1 do if (strmatch(info[i],'*scale*')) then match=i
  if (match ge 0) then begin
    sf = strsplit(info[match],'=',/extract)
    scale = float(sf[1])
    got_scale=1b
  endif else begin
    print,'Base scale not found in dwel header'
    scale=1.0
  endelse
endelse

if (~got_scale) then begin
info_text=[$
'The base fix scale not found in headers',$
'File may be an old version or not correct',$
'Try a new input file ? (Yes/No)']
  result=dialog_message(info_text,/question,$
    title='Casing Maximum not found in headers')
  envi_file_mng,id=fid,/remove
  if (result eq 'Yes') then begin
    goto, infile_select
  endif else begin
    goto, cleanup
  endelse
endif

; Check if file is apparent reflectance
if (DWEL_headers.apprefl_present) then begin
  info_text=[ $
  'File header indicates DWEL data are apparent reflectance data.',$
  'Apparent Reflectance data have already been calibrated and',$
  'attempting to use these data for calibration will not work.',$
  '',$
  'Really use this file for tree calibration ? (No/Yes) ' $
  ]
  result=dialog_message(info_text,/question,$
  title='File is already calibrated!')
  if(result eq 'No') then goto, infile_select
endif

;Get date and time of the acquisition
 match = -1
 for i=0,n_elements(DWEL_headers.DWEL_scan_info)-1 do begin
   if (strmatch(DWEL_headers.DWEL_scan_info[i],'*Data End Time*')) then match=i
 endfor
 if (match ge 0) then begin
   sf = strsplit(DWEL_headers.DWEL_scan_info[match],'=',/extract)
   Data_Time = strtrim(sf[1],2)
 endif else begin
   Data_Time = ''
 endelse

;Get the site description
 match = -1
 for i=0,n_elements(DWEL_headers.DWEL_scan_info)-1 do begin
   if (strmatch(DWEL_headers.DWEL_scan_info[i],'*Scan Description*')) then match=i
 endfor
 if (match ge 0) then begin
   sf = strsplit(DWEL_headers.DWEL_scan_info[match],'=',/extract)
   if (n_elements(sf) gt 1) then begin
     Description_record = strtrim(sf[1],2)
   endif else begin
     Description_record = ''
   endelse
 endif else begin
   Description_record = ''
 endelse
 
 out_scale=1.0

info=DWEL_headers.dwel_adaptation
;now get the DWEL wavelength
match = -1
for i=0,n_elements(info)-1 do BEGIN
   if (strmatch(info[i],'*Wavelength=*', /fold_case)) then match=i
ENDFOR 
IF match GE 0 THEN BEGIN
   text=strtrim(info[match],2)
   kpos=strpos(text,'=')
   wavelength=fix(strtrim(strmid(text,kpos+1,4),2))
ENDIF ELSE BEGIN
    print,strtrim('Input file does NOT have a DWEL file wavelength',2)
    print,'File='+strtrim(infile,2)
    envi_file_mng,id=fid,/remove
    err_flag=1b
    goto,cleanup
ENDELSE 

;check if file has been projected
if (~DWEL_headers.proj_present) then begin
     projected=0b
     Projection_Type='None'
endif else begin
     projected=1b
 match = -1
 for i=0,n_elements(DWEL_headers.DWEL_projection_info)-1 do begin
   if (strmatch(DWEL_headers.DWEL_projection_info[i],'*type*')) then match=i
 endfor
 if (match ge 0) then begin
   sf = strsplit(DWEL_headers.DWEL_projection_info[match],'=',/extract)
   Projection_Type = strtrim(sf[1],2)
 endif else begin
   Projection_Type='None'
 endelse
endelse

print,'projected=',projected
print,'projection type=',projection_type

mask=bytarr(nsamples,nlines)
zenith=fltarr(nsamples,nlines)
azimuth=fltarr(nsamples,nlines)

if (~projected) then begin
;Case for data which have NOT been projected
;find the ancillary file for mask information

  if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
    anc_name=infile+'_ancillary'
  endif else begin
    anc_name=strmid(infile,0,n_dot)+'_ancillary.img'
  endelse

  if(~file_test(anc_name)) then begin
    message_text=[ $
    'Default ancillary file is not present',$
    'default name: '+strtrim(anc_name,2),$
    'Find an ancillary File ? (Yes) or Re-start (No)?' $
    ]
      result=dialog_message(message_text,/question,$
             title='Default ancillary file is NOT present')
      if(result eq 'No') then goto,infile_select
get_anc:
      file = dialog_pickfile(title='Select ancillary DWEL_file', $
             file=anc_name,path=f_path, /must_exist)
;check for error or cancel button hit
      if (file eq '') then begin
        result=dialog_message('Try again ? (No/Yes) ',/question,$
               title='No file selected or operation cancelled !')
        if(result eq 'No') then goto, cleanup else goto, get_anc
      endif
      anc_name=file
  endif

  envi_open_file, anc_name, r_fid=anc_fid,/no_realize
;check if operation cancelled
  if (anc_fid eq -1) then begin
      result=dialog_message('Try again ? (No/Yes) ',/question,$
      title='Error Opening Ancillary Data File !')
      if(result eq 'No') then goto, cleanup else goto, get_anc
  endif

  envi_file_query, anc_fid, nb=nb_anc, nl=nl_anc, ns=ns_anc

  if ((nl_anc ne nlines) or (ns_anc ne nsamples) or (nb_anc lt 3)) then begin
    result=dialog_message('Try again ? (No/Yes) ',/question,$
      title='Ancillary Data File does NOT conform with current DWEL Cube !')
      envi_file_mng,id=anc_fid,/remove
      if(result eq 'No') then goto, cleanup else goto, infile_select
  endif

; Get the Mask, Zenith, Azimuth from the ancillary file
  dims = [-1, 0, nsamples-1, 0, nlines-1]
  Mask = byte(envi_get_data(fid=anc_fid, dims=dims, pos=6))
  Zenith = float(envi_get_data(fid=anc_fid, dims=dims, pos=7))/10.0
  Azimuth = float(envi_get_data(fid=anc_fid, dims=dims, pos=8))/10.0

; Get the Beam Divergence from the Base_fix info.
  info = DWEL_headers.DWEL_scan_info
  match = -1
  for i=0,n_elements(info)-1 do if (strmatch(info[i],'*Beam Divergence*')) then match=i
  if (match ge 0) then begin
    sf = strsplit(info[match],'=',/extract)
    sf = strsplit(sf[1],'mrad',/extract)
    beam_divergence = float(sf)
  endif else begin
    print,'beam divergence not found in dwel header'
    beam_divergence=5.0
  endelse

print,'beam_divergence (mrad)=',beam_divergence

  anc_pos=[6,7,8]
  pos_mask=where(Mask,num_mask)
  if (num_mask le 0) then begin
    info_text=[strtrim(string('The Mask band is all zero !'),2),$
      strtrim('File Name ='+strtrim(anc_name,2),2)]
    result=dialog_message(info_text,/error,title='Mask band is all Masked!')
    print,'The Mask band is all zero!'
    goto,cleanup
  endif else begin
    print,'Number of valid Mask values ='+strtrim(string(100.0*float(num_mask)/float(nsamples*nlines),$
      format='(f10.2)')+' %',2)
    print,'min zenith=',min(zenith[pos_mask])
    print,'max zenith=',max(zenith[pos_mask])
    print,'min azimuth=',min(azimuth[pos_mask])
    print,'max azimuth=',max(azimuth[pos_mask])
    max_zenith=max(zenith[pos_mask])
  endelse
endif else begin

;now for the projected case
;first find the extrainfo file for mask zenith and azimuth information
  if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
      anc_name=infile+'_extrainfo'
  endif else begin
      anc_name=strmid(infile,0,n_dot)+'_extrainfo.img'
  endelse

  if(~file_test(anc_name)) then begin
      message_text=[ $
        'Default ExtraInfo file is not present',$
        'default name: '+strtrim(anc_name,2),$
        'Find an ExtraInfo File ? (Yes) or Re-start (No)?' $
        ]
      result=dialog_message(message_text,/question,$
        title='Default ExtraInfo file is NOT present')
      if(result eq 'No') then goto,infile_select
    get_extra:
      file = dialog_pickfile(title='Select projected ExtraInfo file', $
            file=anc_name,path=f_path, /must_exist)
    ;check for error or cancel button hit
      if (file eq '') then begin
        result=dialog_message('Try again ? (No/Yes) ',/question,$
        title='No file selected or operation cancelled !')
        if(result eq 'No') then goto, cleanup else goto, get_extra
      endif
      anc_name=file
  endif

  envi_open_file, anc_name, r_fid=anc_fid,/no_realize
  ;check if operation cancelled
  if (anc_fid eq -1) then begin
      result=dialog_message('Try again ? (No/Yes) ',/question,$
      title='Error Opening Ancillary Data File !')
      if(result eq 'No') then goto, cleanup else goto, get_extra
  endif

  envi_file_query, anc_fid, nb=nb_anc, nl=nl_anc, ns=ns_anc

  if ((nl_anc ne nlines) or (ns_anc ne nsamples) or (nb_anc lt 4)) then begin
    result=dialog_message('Try again ? (No/Yes) ',/question,$
      title='ExtraInfo File does NOT conform with current Projected DWEL File !')
      envi_file_mng,id=anc_fid,/remove
      if(result eq 'No') then goto, cleanup else goto, infile_select
  endif

  ; Get the mask, zenith and azimuth bands
  dims = [-1, 0, nsamples-1, 0, nlines-1]
  mask = byte(envi_get_data(fid=anc_fid, dims=dims, pos=3))
  zenith = float(envi_get_data(fid=anc_fid, dims=dims, pos=1))/10.0
  azimuth = float(envi_get_data(fid=anc_fid, dims=dims, pos=2))/10.0

; Get the Beam Divergence from the projection info.
  info = DWEL_headers.DWEL_projection_info
  match = -1
  for i=0,n_elements(info)-1 do if (strmatch(info[i],'*output_resolution_(mrad)*')) then match=i
  if (match ge 0) then begin
    sf = strsplit(info[match],'=',/extract)
    beam_divergence = float(sf[1])
  endif else begin
    print,'output resolution not found in dwel header'
    beam_divergence=5.0
  endelse

print,'beam_divergence (mrad)=',beam_divergence

out_scale=1.0

; Get the Beam Divergence from the projection info.
  info = DWEL_headers.DWEL_projection_info
  match = -1
  for i=0,n_elements(info)-1 do if (strmatch(info[i],'*Output scale*')) then match=i
  if (match ge 0) then begin
    sf = strsplit(info[match],'=',/extract)
    out_scale = float(sf[1])
  endif else begin
    print,'output scale from projection not found in dwel header'
    out_scale=1.0
  endelse

print,'beam_divergence (mrad)=',beam_divergence

  anc_pos=[3,1,2]

  pos_mask=where(Mask,num_mask)
  if (num_mask le 0) then begin
    info_text=[strtrim(string('The Mask band is all zero !'),2),$
      strtrim('File Name ='+strtrim(anc_name,2),2)]
    result=dialog_message(info_text,/error,title='Mask band is all Masked!')
    print,'The Mask band is all zero!'
    goto,cleanup
  endif else begin
    print,'Number of valid Mask values ='+strtrim(string(100.0*float(num_mask)/float(nsamples*nlines),$
      format='(f10.2)')+' %',2)
    print,'min zenith=',min(zenith[pos_mask])
    print,'max zenith=',max(zenith[pos_mask])
    print,'min azimuth=',min(azimuth[pos_mask])
    print,'max azimuth=',max(azimuth[pos_mask])
    max_zenith=max(zenith[pos_mask])
  endelse
endelse
zenith=0b
azimuth=0b
mask=0b

;==============================================================================
;set up some defaults and then go and get run parameters
save_rois=0
calibrate=0
fwhm_threshold=2.5
mean_tree_DBH=0.3
tree_refl=0.4

;set up a structure and push it onto the heap
sav={ $
     fwhm_threshold:fwhm_threshold,$
     mean_tree_DBH:mean_tree_DBH,$
     F:F,$
     tree_refl:tree_refl,$
     save_rois:save_rois,$
     calibrate:calibrate,$
     top:event.top $
     }
;now locate the data on the heap with a pointer
p_stat=ptr_new(sav,/no_copy)

;call to widget menu function
status=set_tree_cal_params(p_stat)

if (~status) then begin
  ptr_free, p_stat
  result=dialog_message('No Parameters set in Tree Calibrate - Restart ? ',$
  /question, title='Parameter setup cancelled', $
  dialog_parent=event.top)
  if(result eq 'No') then begin
    goto, cleanup
  endif else goto, input
endif

;set the (valid) sizes back into the default
fwhm_threshold=(*p_stat).fwhm_threshold
mean_tree_DBH=(*p_stat).mean_tree_DBH
save_rois=(*p_stat).save_rois
F=(*p_stat).F
tree_refl=(*p_stat).tree_refl
calibrate=(*p_stat).calibrate

ptr_free,p_stat

; **************
;get the ROI file
if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
  roi_name=infile+'.roi'
endif else begin
  roi_name=strmid(infile,0,n_dot)+'.roi'
endelse

roi_select:
roi_file=envi_pickfile(default=roi_name, filter='*.roi', $
  title='Select DWEL ROI file tree calibration')
if (roi_file eq '') then begin
  result=dialog_message('Try again ? (No/Yes) ',/question,$
  title='No ROI file selected or operation cancelled !')
  if(result eq 'No') then goto, cleanup else goto, roi_select
endif

print,'roi_name=',strtrim(roi_name,2)
print,'roi_file=',strtrim(roi_file,2)

envi_restore_rois,roi_file
cal_roi_ids=envi_get_roi_ids(fid=fid,roi_colors=roi_colors,roi_names=roi_names)
if (cal_roi_ids[0] eq -1) then begin
  print,'cal_roi_ids=',cal_roi_ids
  print,'roi not set up, exit!'
  goto, cleanup
endif
num_rois=n_elements(cal_roi_ids)
print,'number of ROIs=',num_rois
print,'id,color,name'
for i=0,num_rois-1 do begin
  print,cal_roi_ids[i],'  ',roi_colors[i],'  ',strtrim(roi_names[i],2)
endfor

help,range

rstep=abs(range[1]-range[0])

; **************

;Now set up the log file and write out the information obtained so far
; Get path and file name as separate strings

last=strpos(infile,path_sep(),/reverse_search)
in_path=file_dirname(infile)
in_base=strtrim(strmid(infile,last+1,strlen(infile)-last-1),2)

last=strpos(anc_name,path_sep(),/reverse_search)
anc_path = file_dirname(anc_name)
anc_base=strtrim(strmid(anc_name,last+1,strlen(anc_name)-last-1),2)
;
if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
  log_file=infile+'_tree_cal_info.log'
endif else begin
  log_file=strmid(infile,0,n_dot)+'_tree_cal_info.log'
endelse
;see if the log file exists & remove if it does!
if(file_test(log_file)) then begin
  fids=envi_get_file_ids()
  if(fids[0] eq -1) then begin
    file_delete, log_file
  endif else begin
    for i=0,n_elements(fids)-1 do begin
      envi_file_query,fids[i],fname=tname
      if (strlowcase(strtrim(log_file,2)) eq $
          strlowcase(strtrim(tname,2))) then begin
          envi_file_mng,id=fids[i],/remove
      endif
    endfor
    file_delete, log_file
  endelse
endif
;Open Log file
text_err=0
openw, tfile, log_file,/get_lun,error=text_err
if (text_err ne 0) then begin
  print,' '
  print,'error opening the log file'
  print,'Logfile name='+strtrim(log_file,2)
  print,'DWEL_import_batch terminating'
  print,' '
  goto, cleanup
endif
;
printf,tfile,strtrim('DWEL tree calibration Log File',2)
printf,tfile,strtrim('Run made at: '+systime(),2)
printf,tfile,'Basic Run Information'
flush,tfile

;Put info into the log file
printf,tfile,'Description='+strtrim(string(description_record),2)
printf,tfile,'Input_Path='+strtrim(string(in_path),2)
printf,tfile,'Input_File='+strtrim(string(in_base),2)
printf,tfile,'Samples='+strtrim(string(nsamples),2)
printf,tfile,'Lines='+strtrim(string(nlines),2)
printf,tfile,'Ranges='+strtrim(string(nbands),2)
printf,tfile,'Acquisition_Date_Time='+strtrim(string(Data_Time),2)
printf,tfile,'Base_Scale='+strtrim(string(scale),2)
printf,tfile,'Output_Scale='+strtrim(string(out_scale),2)
printf,tfile,'Ancillary_Path='+strtrim(string(anc_path),2)
printf,tfile,'Ancillary_File='+strtrim(string(anc_base),2)
printf,tfile,'Projection='+strtrim(string(Projection_Type),2)
printf,tfile,'Max Range(m)='+strtrim(string(range[nbands-1],format='(f10.2)'),2)
printf,tfile,'Range_Step(m)='+strtrim(string(rstep,format='(f10.4)'),2)
printf,tfile,'Normal Tree Reflectance='+strtrim(string(tree_refl,format='(f10.3)'),2)
printf,tfile,'FWHM Threshold(m)='+strtrim(string(fwhm_threshold,format='(f10.2)'),2)
printf,tfile,'Mean Tree DBH (m)='+strtrim(string(mean_tree_DBH,format='(f10.3)'),2)
printf,tfile,'Beam Divergence (mrad)='+strtrim(string(beam_divergence,format='(f10.3)'),2)
printf,tfile,'Casing Maximum='+strtrim(string(CM,format='(f10.3)'),2)
printf,tfile,'Laser Cal Factor(F)'+strtrim(string(F,format='(f10.4)'),2)
printf,tfile,'calibrate='+strtrim(string(calibrate),2)
flush,tfile

width=fltarr(num_rois)
sinth=fltarr(num_rois)

printf,tfile,'ROI geometric information'
printf,tfile,'ROI_Name,Mean_Mask,Mean_Zenith(deg),Sine_Th,Mean_Azimuth(deg),ROI_Width(mrad)'
print,'ROI_Name,Mean_Mask,Mean_Zenith(deg),Sine_Th,Mean_Azimuth(deg),ROI_Width(mrad)'
for i=0,num_rois-1 do begin
  result=float(envi_get_roi_data(cal_roi_ids[i],fid=anc_fid,pos=anc_pos))
  nsamp=(size(result))[2]
  mean_roi=total(result,2)/float(nsamp)
  min_roi=min(result,dimension=2,max=max_roi)

  width[i]=beam_divergence+100.0*!dtor*(max_roi[2]-min_roi[2])
  sinth[i]=sin(!dtor*mean_roi[1]/10.0)
;
  outstring=strtrim(roi_names[i],2)+','+$
      strtrim(string(mean_roi[0],format='(f10.2)'),2)+','+$
      strtrim(string(mean_roi[1]/10.0,format='(f10.2)'),2)+','+$
      strtrim(string(sinth[i],format='(f10.4)'),2)+','+$
      strtrim(string(mean_roi[2]/10.0,format='(f10.2)'),2)+','+$
      strtrim(string(width[i],format='(f10.3)'),2)
  print,strtrim(outstring,2)
  printf,tfile,strtrim(outstring,2)
endfor
flush,tfile
result=0b

;==============================================
;Now if desired open up the spectral library file
;for the means after removing bad FWHM data

if (save_rois lt 1) then begin
  if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
    outfile=infile+'_tree_cal_means.sli'
  endif else begin
    outfile=strmid(infile,0,n_dot)+'_tree_cal_means.sli'
  endelse

;see if the output file exists and delete if it does
  if(file_test(outfile)) then begin
    fids=envi_get_file_ids()
    if(fids[0] eq -1) then begin
      file_delete, outfile
    endif else begin
      for i=0,n_elements(fids)-1 do begin
        envi_file_query,fids[i],fname=tname
        if (strtrim(strlowcase(outfile),2) eq $
            strtrim(strlowcase(tname),2)) then begin
            envi_file_mng,id=fids[i],/remove
        endif
      endfor
      file_delete, outfile
    endelse
  endif

;Open file for Spectral Library
  text_err=0
  openw, ofile, outfile,/get_lun,error=text_err
  if (text_err ne 0) then begin
    info_text=[ $
    'Error opening output spectral library file '+strtrim(outfile,2),$
    ' ',$
    'Error Name is: '+strtrim(!error_state.name,2)]
    result=dialog_message(info_text, $
     /error,title='Halting file output')
    goto, cleanup
  endif
endif

;==============================================
;Here goes getting information needed from the file about DWEL returns
;then (if set to calibrate) to calculate the calibrations

printf,tfile,'Sample FWHM Values (metres) in ROIs'
printf,tfile,'ROI Name, FWHM values (1 for each shot in ROI)
all_pos=findgen(nbands)
pos_ret=where(range gt 0.25,nranges)

for i=0,num_rois-1 do begin
;get data from input file and samples in the ROI
  result=float(envi_get_roi_data(cal_roi_ids[i],fid=fid,pos=all_pos))
;
;Size[1] is the number of points in an dwel shot
;Size[2] is the number of sample shots in the ROI
;So result is fltarr(nbands,nsamp)
;
  nsamp=(size(result))[2]
  fwhm=rstep*total(result[pos_ret,0:nsamp-1],1)/max(result[pos_ret,0:nsamp-1],dimension=1)
  join=strtrim(fwhm,2)
  join=strjoin(join,',')
  printf,tfile,strtrim(roi_names[i],2)+','+strtrim(join,2),format='(a)'
;deallocate result
  result=0b
  fwhm=0b
  join=0b
endfor
flush,tfile

print,'ROI Name,Nsamp,N_used,FWHM,range,value,F,k(r),pg_r,r_crit'
printf,tfile,'Information for calibration model'
printf,tfile,'ROI Name,Nsamp,N_used,FWHM,range,value,F,k(r),pg_r,r_crit'
flush,tfile

;now get the mean ROI and other statistics
n_rois=0
spec_nam=strarr(num_rois)
xvec=fltarr(num_rois)
yvec=fltarr(num_rois)
for i=0,num_rois-1 do begin
;get data from input file and samples in the ROI
  result=float(envi_get_roi_data(cal_roi_ids[i],fid=fid,pos=all_pos))
;
;Size[1] is the number of points in dwel shot
;Size[2] is the number of samples in the ROI
;So result is fltarr(nbands,nsamp)
;
  nsamp=(size(result))[2]
  fwhm=rstep*total(result[pos_ret,0:nsamp-1],1)/max(result[pos_ret,0:nsamp-1],dimension=1)
  pos_f=where(fwhm lt fwhm_threshold,n_pos_f)
  if (n_pos_f le 0) then begin
    print,'no samples in ROI '+strtrim(string(i+1),2)+' have a valid FWHM'
    goto,move_on
  endif
  n_rois=n_rois+1
  spec_nam[n_rois-1]=strtrim(roi_names[i],2)
  result=result[pos_ret,*]
  result=reform(result[*,pos_f])
  if (n_pos_f gt 1) then begin
    tree_means=total(result,2)/float(n_pos_f)
  endif else tree_means=result
;
  if (save_rois lt 1) then writeu,ofile,tree_means
  max_val=max(tree_means,max_pos)
  r_max=(reform(range[pos_ret]))[max_pos]
  m_fwhm=rstep*total(tree_means)/max_val
  k_r=DWEL_eff_nu(wavelength,beam_divergence,r_max)
  phase_t=phase_trunk(width[i],mean_tree_dbh,r_max,r_crit)
  pg_r=sinth[i]*phase_t
  xvec[n_rois-1]=alog(r_max)
  yvec[n_rois-1]=alog((max_val*F)/(k_r*pg_r))
  outstring=strtrim(roi_names[i],2)+','+$
      strtrim(string(nsamp,format='(i10)'),2)+','+$
      strtrim(string(n_pos_f,format='(i10)'),2)+','+$
      strtrim(string(m_fwhm,format='(f10.4)'),2)+','+$
      strtrim(string(r_max,format='(f10.3)'),2)+','+$
      strtrim(string(max_val,format='(f10.3)'),2)+','+$
      strtrim(string(F,format='(f10.4)'),2)+','+$
      strtrim(string(k_r,format='(f10.3)'),2)+','+$
      strtrim(string(pg_r,format='(f10.3)'),2)+','+$
      strtrim(string(r_crit,format='(f10.3)'),2)
  print,strtrim(outstring,2)
  printf,tfile,strtrim(outstring,2)

;deallocate space
  result=0b
  fwhm=0b
  tree_means=0b
move_on:
endfor
flush,tfile

print,'number of useful ROIs=',n_rois
printf,tfile,'Number of useful ROIs='+strtrim(string(n_rois),2)
flush,tfile

if (n_rois le 0) then begin
  info_text=[$
  'There are NO useful (trunk hit) ROIs in these data!',$
  'Exiting Tree_Cal']
  result=dialog_message(info_text,/error,title='No Useful ROIs!')
  goto,cleanup
endif
spec_nam=spec_nam[0:n_rois-1]
xvec=xvec[0:n_rois-1]
yvec=yvec[0:n_rois-1]
;now if rois being saved set up header etc
if (save_rois lt 1) then begin
  free_lun,ofile,/force
  descrip='Tree Cal Mean Spectra for '+strtrim(infile,2)
  bnames=['Tree Cal Mean Spectra']
  envi_setup_head,fname=outfile,$
    ns=nranges, nl=n_rois,nb=1, file_type=4, interleave=0, $
    data_type=4, wl=range[pos_ret], bnames=bnames, $
    spec_names=spec_nam,descrip=descrip, $
    zplot_titles=['Range(m)','Value'],$
    /write, /open, r_fid=lib_fid
;write out the previous header records
  status=put_headers(lib_fid,DWEL_headers)
endif
;=====================================

;here get the regression and record NSR etc to evaluate fit
;will be checked later in the spreadsheet but first get some parameters

y_bar=mean(yvec)
model=fltarr(n_rois)
result=linfit(xvec,yvec,sigma=sigma,chisq=chisq,yfit=model,prob=prob)

m_bar=mean(model)
Rsq=total((model-replicate(m_bar,n_rois))^2)/ $
total((yvec-replicate(y_bar,n_rois))^2)
print,'RSQ for Model=',Rsq
a=result[0]
cv_a=100.0*sigma[0]/abs(a)
b=result[1]
cv_b=100.0*sigma[1]/abs(b)
RMS=sqrt(total((yvec-model)^2)/float(n_rois))
Cal=exp(a)/(scale*out_scale*tree_refl)
printf,tfile,'Regression for Calibration Fit'
printf,tfile,'a='+strtrim(string(a),2)
printf,tfile,'b='+strtrim(string(b),2)
printf,tfile,'RSQ='+strtrim(string(Rsq),2)
printf,tfile,'RMS='+strtrim(string(RMS),2)
printf,tfile,'CV_a(%)='+strtrim(string(cv_a),2)
printf,tfile,'CV_b(%)='+strtrim(string(cv_b),2)
printf,tfile,'Cal='+strtrim(string(Cal),2)
flush,tfile
printf,tfile,'Regression Data'
printf,tfile,'Number,Xvec,Yvec,Model,Error'
for i=0,n_rois-1 do begin
outstring=strtrim(strtrim(string(i+1),2)+','+ $
  strtrim(string(xvec[i]),2)+','+ $
  strtrim(string(yvec[i]),2)+','+ $
  strtrim(string(model[i]),2)+','+ $
  strtrim(string(yvec[i]-model[i]),2),2)
  printf,tfile,outstring
endfor
flush,tfile
free_lun,tfile,/force
print,'DWEL_tree_cal Finished'
; **************
cleanup:
help,/memory
;clean up any rois
roi_ids = envi_get_roi_ids()
if (roi_ids[0] ge 0)then begin
  print,'number of rois removed=',n_elements(roi_ids)
  envi_delete_rois,/all
endif
;free units and pointers
free_lun,tfile,/force
free_lun,ofile,/force
if(ptr_valid(pb_stats)) then ptr_free,pb_stats
if(ptr_valid(p_stat)) then ptr_free,p_stat
;get state pointer
widget_control,event.top,get_uvalue=pstate
;clean up pointers
widget_control,event.top,/destroy
heap_gc,/verbose
return
;
end
;======================================================================

function DWEL_tree_cal_bannertext
compile_opt idl2

;Provides the banner text as a string array

text=[ $
'Module Description:',$
' ',$
'Procedure DWEL_tree_cal develops calibration factors for DWEL based',$
'on samples of tree returns at varying ranges from the instrument.',$
'The requirements on the target samples are included in a separate',$
'document. The sample locations and extent are indicated by a set',$
'ENVI Regions of Interest (ROIs). The Tree_Cal program takes these',$
'ROIs and a spatially conforming image and extracts the waveforms',$
'in each ROI. They are checked against FWHM and those that are hard',$
'hits are averaged.',$
'',$
'The waveforms should have a single pulse return. The range to the',$
'peak and peak value are extracted. Geometric factors due to the',$
'zenith angle and trunk curvature are taken into account and combined',$
'with the estimated trunk normal reflectance to provide calibrations',$
'to apparent reflectance. Although any heights of trunks may be used',$
'it is best if the ROIs are narrow, not too tall and near the centre',$
'of the tree and near the horizontal (zenith angle near 90 degrees).',$
'',$
'It is usual to run the program first with ROI mean waveforms being',$
'saved but calibration not being attenpted. The data can be analysed',$
'and checked using the output files. When all is well, the program',$
'can be run with the option to calibrate. These results are saved in',$
'a comma separated value file that imports easily into a spreadsheet.',$
'',$
'The program can be run with data after base fixing and saturation',$
'fixing. However, the ROIs must conform with the file used. It is',$
'normal but not necessary to use filtered and projected data. The',$
'program can also be pushed to run on apparent reflectance data to test',$
'the results.',$
'' $
]

return,text

end

pro DWEL_tree_cal_action_go, event
compile_opt idl2

;Event handler for hitting the "GO" button
;on the top level main widget
;simply calls the main action routine

widget_control,event.top,get_uvalue=pstate

(*pstate).go=1

DWEL_tree_cal_doit, event

return

end

pro DWEL_tree_cal_action_exit, event
compile_opt idl2

;Event handler for hitting the "EXIT" button
;on the top level main widget
;cleans up and destroys the widget hierarchy

widget_control,event.top,get_uvalue=pstate
ptr_free, pstate
widget_control,event.top,/destroy

;some protective programming for now
;when all has been well for a while it
;can go away
heap_gc,/verbose

return

end

pro DWEL_tree_cal_resize, event
compile_opt idl2

;Window resize event handler
;for the top level main widget

;First setup protective aunty-Catch to watch out for errors
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
   result=dialog_message(info_text,/error,title='Error in DWEL_tree_cal_resize', $
   dialog_parent=event.top)
   goto, cleanup
endif

;get the information from event - a structure and
;the uvalue - a pointer

 widget_control,event.top,get_uvalue=pstate
 widget_control,event.top, xsize=event.x,ysize=event.y

 xs=max([53*(event.x/(*pstate).xy_old[0]),1])
 ys1=max([10*(event.y/(*pstate).xy_old[1]),1])

 if((event.y/(*pstate).xy_old[1]) eq 1) then begin
  ys2=1
 endif else begin
  ys2=max([2*(event.y/(*pstate).xy_old[1]),1])
 endelse

 xs_act=max([(*pstate).xy_act[0]*(event.x/(*pstate).xy_old[0]),1])
 ys_act=min([max([(*pstate).xy_act[1]*(event.y/(*pstate).xy_old[1]),1]),35])

 widget_control,(*pstate).w_2,xsize=xs,ysize=ys1
 widget_control,(*pstate).w_3,xsize=xs,ysize=ys2
 widget_control,(*pstate).wb_action,xsize=xs_act,ysize=ys_act

 widget_control,event.top,/realize

cleanup:

return

end

pro DWEL_tree_cal_cleanup,top
compile_opt idl2

;The cleanup routine for Xmanager
;This is part of the IDL style top widget

;get state pointer
widget_control,top,get_uvalue=pstate

;clean up pointers
ptr_free, pstate

heap_gc,/verbose

return

end

pro DWEL_tree_cal, event
compile_opt idl2

;The main routine that sets up the widget hierarchy and
;manages the events
;There are only two main events here so it is not so hard
;Providing the overall procedure description in a slider
;menu and the disclaimer are the main purposes
;it is a general top level routine and uses IDL widgets
;for the most part and no auto-manage

; Setup aunty-Catch to watch out for nasty errors
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
   result=dialog_message(info_text,/error,title='Error in DWEL_tree_cal')
   goto, cleanup
endif

;now set up the top level base widget
envi_center,xoff,yoff
wb_base=widget_base(title='Extract Tree Calibration from DWEL files',/column,$
xoffset=xoff,yoffset=yoff,/frame,$
/tlb_size_events)

rev=2s

;get the procedure description and disclaimer to display
intro_text=DWEL_tree_cal_bannertext()
dis_claim=DWEL_disclaimer(rev)

;Set up a base widget to hold the description and disclaimer
wb_info=widget_base(wb_base,/column)

;Use two slider label box routines from the ENVI library
w_2=widget_slabel(wb_info,prompt=intro_text,/frame,xsize=53,ysize=10)
w_3=widget_slabel(wb_info,prompt=dis_claim,/frame,xsize=53,ysize=1)

;Now set up another base widget for the action buttons
wb_action=widget_base(wb_base,/row)

;Action buttons are just "go" and "exit"
w_go=widget_button(wb_action,value='Go',uvalue='DWEL_tree_cal_action_go',$
event_pro='DWEL_tree_cal_action_go',/dynamic_resize)
w_exit=widget_button(wb_action,Value='Exit',uvalue='DWEL_tree_cal_action_exit',$
event_pro='DWEL_tree_cal_action_exit',/dynamic_resize)

;realise the widget hierarchy
widget_control,wb_base,/realize

go=0

widget_geometry_1=widget_info(wb_base,/geometry)
widget_geometry_2=widget_info(wb_action,/geometry)

xy_old=[widget_geometry_1.xsize,widget_geometry_1.ysize]
xy_act=[widget_geometry_2.xsize,widget_geometry_2.ysize]

;set up a state structure
state={ $
tlb:wb_base,$
w_2:w_2, $
w_3:w_3, $
wb_action:wb_action, $
go:go, $
xy_old:xy_old, $
xy_act:xy_act, $
wb_info:wb_info $
}

pstate=ptr_new(state,/no_copy)
widget_control,wb_base,set_uvalue=pstate

;call xmanager
xmanager,'DWEL_tree_cal',wb_base, $
event_handler='DWEL_tree_cal_resize',$
cleanup='DWEL_tree_cal_cleanup',/no_block

;exit back to ENVI
;leave a bit of a trace for bad exits for now
;this can go when things stabilise - along with
;much protective programming and aunty catch
;in all the many places it is at this time!
cleanup:
heap_gc,/verbose
return

end

