; Simple AT projection of unprocessed data cube.
; Arrange pixels in array of zenith and azimuth angle instead of shot number and scan number

pro dwel_val_block
  pi_val=4.0*atan(1.0)
  rad=pi_val/180.0
  deg=1.0/rad
  eta=1.0e-7
end

; Max_Zenith_Angle: in unit of degree
; output_resolution: in unit of mrad
; overlap: azimuth range of overlapping area, in unit of degree
pro dwel_cube2at_nsf, DWEL_Cube_File, DWEL_Anc_File, DWEL_AT_File, $
  Max_Zenith_Angle, output_resolution, zen_tweak, err, Overlap=overlap

  compile_opt idl2
  envi, /restore_base_save_files
  envi_batch_init, /no_status_window

  print, 'entering dwel_cube2at'
    
  inlun=100
  ofile=30
  p_stat=0
  p_list=0
  err=0
  err_flag=0b
  
  ;;;;;;;;;;;;;;;;;;;;;;;;
  ; some default paramters
  def_ifov_x=4.0
  ifov_min=0.0
  ifov_max=60.0
  tmax_min=0.0
  tmax_max=180.0
  def_tmax=110.0
  projection_index=0
  type=4
  scale=1.0d0
  r2mr=1000.0
  Proj_name=['Hemispherical','Andrieu Normal','Andrieu Transpose']
  ;set speed of light metres per nsec /2
  c=0.299792458
  c2=c/2.0
  ;;;;;;;;;;;;;;;;;;;;;;;;
  
  envi_open_file, DWEL_Cube_File,r_fid=fid,/no_interactive_query,/no_realize
  if(fid eq -1) then begin
    print, 'Processing stopped! Failed to open the input data cube file: ' + DWEL_Cube_File
    err=1
    err_flag=1b
    goto, cleanup
  endif
  
  ;get the input image dimensions and other info
  envi_file_query, fid, ns=ns, nl=nl, nb=nb, $
    byte_swap=order, data_type=type, interleave=ftype, $
    bnames=band_name, dims=dims, file_type=f_type, $
    xstart=xstart, ystart=ystart,wl=wl
    
  result=envi_file_type(f_type)
  
  print,'input file type = ',result
  
  if (result ne 'ENVI Standard') then begin
    print, 'Processing stopped! Input file not ENVI Standard!'
    print,'file_type='+strtrim(result)
    envi_file_mng,id=fid,/remove
    err=2
    err_flag=1b
    goto, cleanup
  endif
  
  samples=ns
  lines=nl
  bands=nb
  
  ;get number of bytes in the input data
  nbytes=dt2nb(type)
  
  band_pos=indgen(nb)
  help,wl
  
  x_range=[dims[1],dims[2]]
  y_range=[dims[3],dims[4]]
  
  ;set the type of file
  ft_nam='Unknown'
  case ftype of
    0: ft_nam='BSQ'
    1: ft_nam='BIL'
    2: ft_nam='BIP'
  endcase
  
  ;set the data type
  dt_nam='Unknown'
  case type of
    1: dt_nam='Byte'
    2: dt_nam='Int'
    3: dt_nam='Long Int'
    4: dt_nam='Floating Point'
    5: dt_nam='DP Floating Point'
    6: dt_nam='Complex'
    9: dt_nam='DP Complex'
    12: dt_nam='Unsigned Int'
    13: dt_nam='Unsigned Long Int'
    14: dt_nam='64-bit Int'
    15: dt_nam='64-bit Long Int'
  endcase
  
  ;get path and image name as separate strings
  last=strpos(DWEL_Cube_File,path_sep(),/reverse_search)
  f_path=strmid(DWEL_Cube_File,0,last+1)
  f_base=strtrim(strmid(DWEL_Cube_File,last+1,strlen(DWEL_Cube_File)-last-1),2)
  
  ;check for the ancillary data file
  n_base=strlen(DWEL_Cube_File)
  n_dot=strpos(DWEL_Cube_File,'.',/reverse_search)
  
  ;now get the DWEL headers that are present
  ;set up a base structure for the DWEL headers
  DWEL_headers={ $
    f_base:f_base $
    }
    
  ;find all of the DWEL headers in the hdr file as defined by FID
  status=DWEL_get_headers(fid,DWEL_headers)
  
  if (not status) then begin
    print, 'Processing stopped! Bad FID in DWEL_get_headers on input file!'
    envi_file_mng,id=fid,/remove
    err=3
    err_flag=1b
    goto, cleanup
  endif
  
  if (DWEL_headers.headers_present le 0s or not DWEL_headers.run_present) then begin
    print,'Processing stopped! Input File is NOT a valid DWEL Cube file!'
    print,'DWEL_headers.headers_present='+strtrim(string(DWEL_headers.headers_present),2)
    print,'DWEL_headers.run_present='+strtrim(string(DWEL_headers.run_present),2)
    envi_file_mng,id=fid,/remove
    err=4
    err_flag=1b
    goto, cleanup
  endif
  
  ;locate and set the scale factor for the average images
  if (not DWEL_headers.base_present) then begin
    base_scale=1.0
    if (type le 1) then begin
      scale=1000.0d0
    endif else begin
      scale=100.0d0
    endelse
  endif else begin
    text=strtrim(DWEL_headers.DWEL_base_fix_info[8],2)
    l=strlen(text)
    k=strpos(text,'=')
    if (strtrim(strmid(text,0,k),2) ne 'scale') then begin
      base_scale=1.0
      if (type le 1) then begin
        scale=1000.0d0
      endif else begin
        scale=100.0d0
      endelse
    endif else begin
      reads,strtrim(strmid(text,k+1,l),2),base_scale
      if (type le 1) then begin
        scale=1000.0d0/double(base_scale)
      endif else begin
        scale=100.0d0/double(base_scale)
      endelse
    endelse
  endelse
  
  ;Get date and time of the acquisition
  DWEL_date_time=''
  match = -1
  for i=0,n_elements(DWEL_headers.DWEL_scan_info)-1 do begin
    if (strmatch(DWEL_headers.DWEL_scan_info[i],'*Data End Time*') or $
      strmatch(DWEL_headers.DWEL_scan_info[i],'*DWEL Date Time*')) then match=i
  endfor
  if (match ge 0) then begin
    sf = strsplit(DWEL_headers.DWEL_scan_info[match],'=',/extract)
    DWEL_date_time = strtrim(strcompress(sf[1]),2)
  endif else begin
    DWEL_date_time = ''
  endelse
  
  ;Get the site description
  DWEL_description_record=''
  match = -1
  for i=0,n_elements(DWEL_headers.DWEL_scan_info)-1 do begin
    if (strmatch(DWEL_headers.DWEL_scan_info[i],'*Scan Description*')) then match=i
  endfor
  if (match ge 0) then begin
    sf = strsplit(DWEL_headers.DWEL_scan_info[match],'=',/extract)
    if (n_elements(sf) gt 1) then begin
      DWEL_description_record = strtrim(sf[1],2)
    endif else begin
      DWEL_description_record = ''
    endelse
  endif else begin
    DWEL_description_record = ''
  endelse
  
  print,''
  print,'DWEL run Date & Time='+strtrim(DWEL_date_time,2)
  print,'DWEL run Description='+strtrim(DWEL_description_record,2)
  
  if (not DWEL_headers.base_present) then begin
    info = DWEL_headers.DWEL_scan_info
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
    if (~srate_set) then print,'sampling rate not set'
    print,'sampling rate=',srate
    tdif=c2/srate
    print,'time step='+strtrim(string(tdif),2),' metres'
    wl=findgen(bands)*tdif
  endif
  
  if(not file_test(DWEL_Anc_File)) then begin
    print,'Processing stopped! Ancillary file not present or not correct!'
    envi_file_mng,id=fid,/remove
    err=5
    err_flag=1b
    goto, cleanup
  endif
  
  text_err=0
  envi_open_file, DWEL_Anc_File,r_fid=anc_fid,/no_interactive_query,/no_realize
  if (anc_fid eq -1) then begin
    print,'Processing stopped! Error opening ancillary data file '+strtrim(DWEL_Anc_File,2)
    envi_file_mng,id=fid,/remove
    err=6
    err_flag=1b
    goto, cleanup
  endif
  
  ;get the input image dimensions and other info
  envi_file_query, anc_fid, ns=Nshots, nl=Nscans, nb=nb_anc
  
  if ((Nshots ne Samples) or (Nscans ne lines)) then begin
    print,'Processing stopped! Dimension of ancillary File '+strtrim(DWEL_Anc_File,2)+' does NOT Conform with input !'
    envi_file_mng,id=fid,/remove
    envi_file_mng,id=anc_fid,/remove
    err=7
    err_flag=1b
    goto, cleanup
  endif
  
  ;now get the DWEL headers that are present for the ancillary file
  ;set up a base structure for the  headers
  DWEL_anc_headers={ $
    f_base:DWEL_Anc_File $
    }
    
  ;find all of the DWEL headers in the hdr file as defined by FID
  status=DWEL_get_headers(anc_fid,DWEL_anc_headers)
  
  if (not status) then begin
    print,'Processing stopped! Bad FID in DWEL_get_headers for Ancillary File'
    envi_file_mng,id=fid,/remove
    envi_file_mng,id=anc_fid,/remove
    err=8
    err_flag=1b
    goto, cleanup
  endif
  
  if (DWEL_anc_headers.headers_present le 0s or not DWEL_anc_headers.run_present) then begin
    print,'DWEL_anc_headers.headers_present='+strtrim(string(DWEL_anc_headers.headers_present),2)
    print,'DWEL_anc_headers.run_present='+strtrim(string(DWEL_anc_headers.run_present),2)
    print,'Processing stopped! Ancillary file NOT a valid DWEL Cube file'
    envi_file_mng,id=fid,/remove
    envi_file_mng,id=anc_fid,/remove
    err=9
    err_flag=1b
    goto, cleanup
  endif
  
  print,'Input information from input files complete'
  
  ;set output band number
  nb_out=n_elements(where(band_pos gt -1))
  
  ;get the mask
  Mask_all=bytarr(Nshots,Nscans)+1b
  dims=[-1,0,Nshots-1,0,Nscans-1]
  Mask_all=envi_get_data(fid=anc_fid,dims=dims,pos=6)
  
  ;get encoder positions!
  ShotZen=fltarr(Nshots,Nscans)
  ShotAZim=fltarr(Nshots,Nscans)
  ShotZen=float(envi_get_data(fid=anc_fid,dims=dims,pos=2))
  ShotAZim=float(envi_get_data(fid=anc_fid,dims=dims,pos=3))
  
  envi_file_mng,id=anc_fid,/remove
  envi_file_mng,id=fid,/remove
  
  ;; ;now get start of data to ensure wrap
  ;; sel=0
  ;; bad=50
  ;; bad_end=10
  ;; RotEnc_max=524288.0d0
  ;; Rot_med=double(median(ShotAzim,dimension=1))
  ;; ;
  ;; Rot_med=reform(Rot_med[bad:Nscans-1])
  ;; pos=where(Rot_med gt 0,npos)
  ;; if (npos gt 0) then Rot_med=reform(Rot_med[pos])
  ;; numval=n_elements(Rot_med)
  
  ;; Rot_med_compl=Rot_med-RotEnc_Max/2.0d0
  ;; posc=where(Rot_med_compl lt 0.0,nposc)
  ;; if (nposc gt 0) then Rot_med_compl[posc]=Rot_med_compl[posc]+RotEnc_Max
  ;; posc=0b
  
  ;; valend=Rot_med[numval-bad_end-1]
  ;; posc=where(Rot_med_compl lt valend,nposc)
  ;; if (nposc gt 0) then begin
  ;;   sel=max([bad,pos[posc[0]]+bad])
  ;; endif else begin
  ;;   sel=bad
  ;; endelse
  ;; posc=0b
  ;; pos=0b
  ;; print,'final sel=',sel

  bad=50 ;; always discard the first 50 scan lines of possible bad rotary
  ;; encoders due to intertial of the rotation or lack of lock-in of rotary
  ;; encoder.
  bad_end=10 ;; always discard the last 10 scan lines just in case of any funky
  ;; error.
  ;; set the mask for overlapping, 0: pixels to be discarded.
  ;; here we are using ENCODER values, NOT actual angular values.
  ;; by default, i.e. no overlap is given, do NOT remove overlap and leave the
  ;; scan as it is collected
  ;; Now calcualte the azimuth angle of each scan line
  tmpazim = fltarr(nscans)
  for i=0,nscans-1 do begin
    tmp = ShotAzim[where(Mask_all[*, i]), i]
    tmpazim[i] = median(tmp)
  endfor
  ;; if the decreasing rotary encoder passes through 0 and 524288 (2*pi), add a
  ;; round of 524288 (2*pi) so that later angular calculation is easier.
  tmpdiff = tmpazim[0:nscans-2] - tmpazim[1:nscans-1]
  tmppos = where(tmpdiff lt -524288.0d0/2.0, tmpcount)
  tmpazim = tmpazim[0:tmppos[0]] + 524288.0d0
  ;; calculate overlap azimuth range in the scan
  az_range = tmpazim[bad+1] - tmpazim[nscans-1-bad_end]
  scan_overlap = (az_range - 524288.0d0/2.0) / 524288.0d0 * 360.0
  if scan_overlap lt 0 then begin
    scan_overlap = 0.0
  endif
  if n_elements(overlap) ne 0 or arg_present(overlap) then begin
    ;; overlap is given
    ;; set mask according to azimuth angles
    ;; we've found that the azimuth encoders of the first few scan lines
    ;; of a scan could be of no change possibly due to the inertial of the
    ;;instrument rotation or lack of lock-in of rotary encoder. Thus here we
    ;;retain the last 180 degrees of scan lines and discard the first few scan
    ;;lines with possibly bad azimuth values.
    firstazim = tmpazim[nscans-1-bad_end] + 524288.0d0/2.0 + float(overlap)/360.0*524288.0d0
    tmppos = where(tmpazim gt firstazim, tmpcount)
    if tmpcount gt 0 then begin
      if tmppos[tmpcount-1] gt bad then begin
        sel = tmppos[tmpcount-1]
      endif else begin
        sel = bad
      endelse
    endif else begin
      sel = bad
    endelse
    Mask_all[*, 0:sel] = 0
    print, 'first scan line in the projection = ', sel+1
    Mask_all[*, (nscans-bad_end):(nscans-1)] = 0
    print, 'last scan line in the projection = ', nscans-bad_end-1
    print, 'Given overlap azimuth = ', overlap
    print, 'Scan''s own overlap azimuth = ', scan_overlap
    if overlap gt scan_overlap then begin
      overlap = scan_overlap
    endif
  endif else begin
    ;; no overlap is given
    ;; do not change the mask and use all valid pixels in the projection
    ;; calculate the overlapping azimuth angle in the scan
    Mask_all[*, 0:bad] = 0
    print, 'first scan line in the projection = ', bad+1
    Mask_all[*, (nscans-bad_end):(nscans-1)] = 0
    print, 'last scan line in the projection = ', nscans-bad_end-1
    print, 'Project the whole scan, overlap azimuth = ', scan_overlap
    overlap = scan_overlap
  endelse
  
  ;===========================================
  ;set up a structure and push it onto the heap
  sav={ $
    Nshots:Nshots,$
    Nscans:Nscans,$
    ShotZen:ShotZen,$
    ShotAzim:ShotAzim $
    }
    
  ;now put the data on the heap with a pointer
  p_stat=ptr_new(sav,/no_copy)
  
  status = dwel_set_theta_phi_nsf(p_stat,zen_tweak)
  
  ;put the results into the local arrays
  ShotZen=(*p_stat).ShotZen
  ShotAzim=(*p_stat).ShotAzim
  ptr_free, p_stat
  
  ;set the mask for angles outside range
  pos=where(Shotzen gt Max_Zenith_Angle,npos)
  print,'npos angles above max=',npos
  
  if (npos gt 0) then Mask_all[pos]=0
  pos=0b
  
  ;; ;now introduce the new mask to ensure minimum wrap
  ;; mask_all[*,0:sel]=0
  ;; mask_all[*,(nscans-bad_end):(nscans-1)]=0
  
  srate_set=0b
  ;set the default sampling rate
  match = -1
  for i=0,n_elements(DWEL_headers.DWEL_scan_info)-1 do begin
    if (strmatch(DWEL_headers.DWEL_scan_info[i],'*Sampling Rate*')) then match=i
  endfor
  if (match ge 0) then begin
    text=strtrim(DWEL_headers.DWEL_scan_info[match],2)
    k=strpos(text,'=')
    l=strpos(text,'smp/ns')
    var2=float(strtrim(strmid(text,k+1,l-k-1),2))
    if (var2 gt 0.0) then begin
      srate=var2
      srate_set=1
    endif else begin
      srate=2.0
      srate_set=0b
    endelse
  endif else begin
    srate=2.0
    srate_set = 0b
  endelse
  
  if (~srate_set) then print,'Sampling rate NOT read from headers!!!'
  print, 'Sampling rate=' + strtrim(string(srate),2) 
  
  ;Get the beam divergence
  buf=''
  match = -1
  for i=0,n_elements(DWEL_headers.DWEL_scan_info)-1 do begin
    if (strmatch(DWEL_headers.DWEL_scan_info[i],'*Beam Divergence*')) then match=i
  endfor
  if (match ge 0) then begin
    ;   print,'match=',DWEL_headers.DWEL_scan_info[match]
    sf = strsplit(DWEL_headers.DWEL_scan_info[match],'=',/extract)
    if (n_elements(sf) gt 1) then begin
      buf = strtrim(strcompress(sf[1]),2)
    endif else begin
      buf = ''
    endelse
  endif else begin
    buf = ''
  endelse
  l=strpos(buf,'mrad')
  if (l gt 0) then buf=strtrim(strmid(buf,0,l-1),2) else buf=''
  ; print,'buf=',buf
  val=float(buf)
  beam_div=val[0]
  
  ;Get the scan_step
  buf=''
  match = -1
  for i=0,n_elements(DWEL_headers.DWEL_scan_info)-1 do begin
    if (strmatch(DWEL_headers.DWEL_scan_info[i],'*Scans per Complete Rotation*')) then match=i
  endfor
  if (match ge 0) then begin
    ;   print,'match=',DWEL_headers.DWEL_scan_info[match]
    sf = strsplit(DWEL_headers.DWEL_scan_info[match],'=',/extract)
    if (n_elements(sf) gt 1) then begin
      buf = strtrim(strcompress(sf[1]),2)
    endif else begin
      buf = ''
    endelse
  endif else begin
    buf = ''
  endelse
  val=float(buf)
  if (val[0] gt 0.0) then scan_step=2000.0*!pi/val[0] else scan_step=0.0 ; unit of scan step: mrad
  sampling_ratio=beam_div/scan_step
  
  ;get the Ifov and Maximum Theta
  set_ifov:
  
  ifov_x=output_resolution
  t_max=Max_Zenith_Angle*!pi/180.0
  ptype='AT'
  pname='Andrieu Transpose'
  
  print,'ifov_x=',ifov_x
  
  ;now convert the data to (x,y) coordinates
  
  print, 'Andrieu_transpose'
  status=t_Andrieu_tp2xy(Nshots*Nscans,ShotZen,ShotAzim,x_proj,y_proj)
  
  ShotZen=0b
  ShotAzim=0b
  
  ;define image geometry - linear size is (2*k_inc+1) cells
  
  ;set output dimensions
  ns_out=fix(r2mr*t_max/ifov_x)+1
  nl_out=fix(r2mr*2.0*!pi/ifov_x)+1
  
  ;set the step size in (x,y) space corresponding to the resolution
  h=(t_max/float(ns_out))*!radeg
  h21=1.1*(h/2.0)
  h22=1.1*((!radeg*scan_step)/r2mr)/2.0
  h2=max([h21,h22])
  print, 'projection step size: ', h2
  
  counter=0L
  Tot_Count=0L
  gotind=0b
  a_ref=[Nshots,Nscans]
  
  ;Now set up the pointer array to the cell data and all is ready
  
  p_list=make_array(nl_out,/ptr)
  num_val=make_array(ns_out,nl_out,/long)
  theta=make_array(ns_out,nl_out,/long)
  phi=make_array(ns_out,nl_out,/long)
  mask=make_array(ns_out,nl_out,/long)+1L
  
  A_theta=make_array(ns_out,/double)
  A_phi=make_array(nl_out,/double)
  pos_ind=0b
  
  ;; set a scale for angular output
  angle_scale = 100.0
  ;now loop over the cells and set up pointers to the original image cells
  ;involved in the output spatial cells
  for i=0, nl_out-1 do begin
    y=float(i)*h
    y_min=y-h2
    y_max=y+h2
    pos_y=where(((y_proj ge y_min) and (y_proj le y_max) and Mask_all),count_y)
    
    if (count_y gt 0) then begin
      for j=0,ns_out-1 do begin
        x=float(ns_out-1-j)*h
        ;get four corners of the pixel
        x_min=x-h2
        x_max=x+h2
        pos_x=where(x_proj[pos_y] ge x_min and x_proj[pos_y] le x_max,count_x)
        if(count_x gt 0) then begin
          num_val[j,i]=long(count_x)
          temp=array_indices(a_ref,reform(pos_y[pos_x]),/dimensions)
          jj=size(temp,/dimensions)
          if (n_elements(jj) ne 2) then begin
            a_add=j
            temp=[temp,a_add]
          endif else begin
            a_add=intarr(jj[1])+j
            temp=[temp,transpose(a_add)]
          endelse
          
          if (gotind) then begin
            pos_ind=[[pos_ind],[temp]] ; pos_ind is n*3 array, n is the number
                                ; of shots in the line i of the projection,
                                ; first and second columns are the pixel
                                ; location in the original data cube, the third
                                ; line is the column j of the projection.
          endif else begin
            pos_ind=[[temp]]
            gotind=1b
          endelse
          counter=counter+1L
          Tot_Count=Tot_Count+long(count_x)
          temp=0b
        endif else begin
          mask[j,i]=0L
        endelse
        status=t_Andrieu_xy2tp(1,x,y,th,ph)
        A_theta[j]=th
        A_phi[i]=ph
        theta[j,i]=round(angle_scale*th)
        phi[j,i]=round(angle_scale*ph)
      endfor
      if (total(reform(num_val[*,i])) gt 0) then begin
        temp=sort(reform(pos_ind[1,*])) ; sort the pos_ind by column, the shot number/zenith
        pos_ind=pos_ind[*,temp]
      endif
    endif else begin
      for j=0,ns_out-1 do begin
        x=float(ns_out-1-j)*h
        num_val[j,i]=0L
        mask[j,i]=0L
        status=t_Andrieu_xy2tp(1,x,y,th,ph)
        A_theta[j]=th
        A_phi[i]=ph
        theta[j,i]=round(angle_scale*th)
        phi[j,i]=round(angle_scale*ph)
      endfor
    ;print,'Hit an empty line at i=',i
    endelse
    p_list[i]=ptr_new(pos_ind,/no_copy)
    pos_ind=0b
    gotind=0b
    temp=0b
  endfor
  
  a_zen=strtrim(string(A_theta,format='(f12.4)'),2)
  a_azm=strtrim(string(A_phi,format='(f12.4)'),2)
  
  DWEL_Andrieu_zenith=strtrim('zen=('+strcompress(strjoin(a_zen,',',/single),/remove_all)+')',2)
  DWEL_Andrieu_azimuth=strtrim('azm=('+strcompress(strjoin(a_azm,',',/single),/remove_all)+')',2)
  
  a_zen=0b
  a_azm=0b
  A_theta=0b
  A_phi=0b
  
  if (type le 1) then begin
    scaler=100.0d0
  endif else begin
    scaler=1.0d0
  endelse
  DWEL_Projection_info=strarr(14)
  DWEL_Projection_info=[ $
    'Program='+'DWEL_Cube2AT_NSF Projection Routine',$
    'Processing_Date_Time='+strtrim(systime(),2),$
    'Projection_type='+ptype,$
    'Projection_name='+pname,$
    'Beam_Divergence_(mrad)='+strtrim(string(beam_div,format='(f14.3)'),2),$
    'Scan_Step_(mrad)='+strtrim(string(scan_step,format='(f14.3)'),2),$
    'Sampling_Ratio='+strtrim(string(sampling_ratio,format='(f14.3)'),2),$
    'output_resolution_(mrad)='+strtrim(string(ifov_x,format='(f10.2)'),2),$
    'max_zenith_angle_(deg)='+strtrim(string(!radeg*t_max,format='(f10.2)'),2),$
    'Zen_tweak_(enc)='+strtrim(string(zen_tweak),2),$
    'Mean_image_scale='+strtrim(string(scale,format='(f10.2)'),2),$
    'Output_scale='+strtrim(string(scaler,format='(f10.2)'),2), $
    'Angular_scale='+strtrim(string(angle_scale,format='(f10.2)'),2), $
    'Overlap_azimuth='+strtrim(string(overlap, format='(f10.3)'), 2) $
    ]
  DWEL_Projection_info=strtrim(DWEL_Projection_info,2)
  
  ;all ready to go ... so get output file name[s]
  output_envi:
  
  text_err=0
  ;open input file
  openr,inlun,DWEL_Cube_File,/get_lun,error=text_err
  if (text_err ne 0) then begin
    print,'Error opening input file '+strtrim(DWEL_Cube_File,2)
    if (ptr_valid(p_list[0])) then ptr_free,p_list
    err=9
    err_flag=1b
    goto,cleanup
  endif
  ;Open output file
  text_err=0
  openw, ofile, DWEL_AT_File,/get_lun,error=text_err
  if (text_err ne 0) then begin
    print,'Error opening output file '+strtrim(DWEL_AT_File,2)
    if (ptr_valid(p_list[0])) then ptr_free,p_list
    err=10
    err_flag=1b
    goto, cleanup
  endif
  
  print,'pre-processing done - projecting the image!'
  
  ;set the output data type
  if (type lt 4 or type gt 9) then begin
    out_type=2
  endif else begin
    out_type=4
  endelse
  
  ;set up the arrays for data and output
  data=make_array(Nshots,nb_out,/double)
  accum=make_array(ns_out,nl_out,/double)
  accum_abs=make_array(ns_out,nl_out,/double)
  temp=make_array(Nshots,nb_out,/double)
  maxwf = make_array(ns_out, nl_out, /double)
  meanwf = make_array(ns_out, nl_out, /double)
  accum_r=make_array(ns_out,nl_out,/double)
  accum_r2=make_array(ns_out,nl_out,/double)
  
  ;set up the range scaling
  
  sc_r=dblarr(nb_out)
  sc_r2=dblarr(nb_out)
  
  pos_nz=where(wl ge 0.0,n_pos)
  step=wl[1]-wl[0]
  
  if (n_pos gt 0) then begin
    sc_r[pos_nz]=double(step)*(dindgen(n_pos)+1.0d)
    sc_r2[pos_nz]=(double(step)*(dindgen(n_pos)+1.0d))^2
  endif
  pos_nz=0b
  num_avg=make_array(ns_out,nl_out,/long)
  
  ;set the pointers
  bufrs=long64(nbytes)*long64(nb)*long64(ns)
  pointsz=long64(0)
  
  ;do the processing over the output tiles
  ;BIL Tile
  
  for k=0L, nl_out-1 do begin
    current=-1
    ;count=0L
    temp=make_array(ns_out,nb_out,/double)
    ;help,*(p_list[k])
    pos_nz=where(num_val[*,k] gt 0,n_pos)
    if (n_pos gt 0) then begin
      pos_ind=*(p_list[k])
      kk=n_elements(reform(pos_ind[0,*]))
      point=0L
      while (point lt kk) do begin
        lin=pos_ind[1,point]
        if (lin ne current) then begin
          data=0s
          ;          data=double(envi_get_tile(tile_id,lin))
          ;          data=double(envi_get_slice(fid=fid, line=lin, /bil))
          pointsz=long64(lin)*long64(bufrs)
          data=read_binary(inlun,data_start=pointsz,data_dims=[ns,nb],data_type=type)
          data=double(data)
          ;count=count+1L
          current=lin
        endif
        ;do something!
        temp[pos_ind[2,point],*]=temp[pos_ind[2,point],*]+ $
          data[pos_ind[0,point],*]
        num_avg[pos_ind[2,point],k]=num_avg[pos_ind[2,point],k]+1L
        point=point+1L
      endwhile
      temp[pos_nz,*]=cmreplicate((1.0/float(reform(num_avg[pos_nz,k]))),nb_out)*temp[pos_nz,*]
    endif
    pos_ind=0b
    pos_nz=0b
    data=0b
    writeu,ofile,fix(round(scaler*temp))
    
    ;now get the statistics of the image for the extra info file
    accum[*,k]=total(abs(temp),2,/double)/float(nb_out)
    accum_abs[*,k]=total(abs(temp),2,/double)/float(nb_out)
    accum_r[*,k]=transpose(abs(temp))##sc_r
    accum_r2[*,k]=transpose(abs(temp))##sc_r2
    
    maxwf[*,k] = max(abs(temp), dimension = 2)
    meanwf[*,k] = mean(temp, dimension = 2)
    temp=0b
  endfor
  
  ;  envi_tile_done, tile_id
  free_lun,inlun,/force
  free_lun,ofile,/force
  ptr_free, p_list
  data=0b
  temp=0b
  
  if (total(abs(num_avg-num_val)) gt 0) then begin
    print,'numbers of averaged cells do NOT agree!'
    print,'total counted in first loop=',total(num_val)
    print,'total counted in second loop=',total(num_avg)
  endif
  
  print,'Done projecting - now get info and clean up!'
  
  ;; compute the stars
  star_r=make_array(ns_out,nl_out,/double)
  star_r2=make_array(ns_out,nl_out,/double)
  pos=where(abs(accum_abs) gt 1.0e-5,count)
  if (count gt 0) then begin
    star_r[pos]=abs(accum_r[pos])/abs(accum_abs[pos])
    star_r2[pos]=sqrt(abs(accum_r2[pos])/abs(accum_abs[pos]))
  endif
  accum_abs=0b
  pos=0b
  
  accum_r=accum_r/total(sc_r)
  accum_r2=accum_r2/total(sc_r2)
  
  sc_r=0b
  sc_r2=0b
  
  image_statistics,maxwf,mask=mask,minimum=amin,maximum=amax,mean=amean,stddev=asdev
  scaler=4095.0/(amax-amin)
  maxwf=((round(scaler*(maxwf-amin))>0L)<4095L)
  
  image_statistics,accum,mask=mask,minimum=amin,maximum=amax,mean=amean,stddev=asdev
  scaler=4095.0/(amax-amin)
  accum=((round(scaler*(accum-amin))>0L)<4095L)
  ;
  ;  print,'accum mean=',amean
  
  DWEL_projection_info=[DWEL_projection_info, $
    'Stats_Format=(Min,Mean,Max,Stddev)', $
    'accum_Stats=('+strtrim(string(amin),2)+',' $
    +strtrim(string(amean),2)+',' $
    +strtrim(string(amax),2)+',' $
    +strtrim(string(asdev),2)+')' $
    ]
  ;
  image_statistics,accum_r,mask=mask,minimum=amin,maximum=amax,mean=amean,stddev=asdev
  scaler=4095.0/(amax-amin)
  accum_r=((round(scaler*(accum_r-amin))>0L)<4095L)
  ;
  DWEL_projection_info=[DWEL_projection_info, $
    'accum_r_Stats=('+strtrim(string(amin),2)+',' $
    +strtrim(string(amean),2)+',' $
    +strtrim(string(amax),2)+',' $
    +strtrim(string(asdev),2)+')' $
    ]
  ;
  image_statistics,accum_r2,mask=mask,minimum=amin,maximum=amax,mean=amean,stddev=asdev
  scaler=4095.0/(amax-amin)
  accum_r2=((round(scaler*(accum_r2-amin))>0L)<4095L)
  ;
  DWEL_projection_info=[DWEL_projection_info, $
    'accum_r2_Stats=('+strtrim(string(amin),2)+',' $
    +strtrim(string(amean),2)+',' $
    +strtrim(string(amax),2)+',' $
    +strtrim(string(asdev),2)+')' $
    ]
  ;
  ;  image_statistics,star_r,mask=mask,minimum=amin,maximum=amax,mean=amean,stddev=asdev
  scaler=4095.0/(amax-amin)
  star_r=((round(scaler*(star_r-amin))>0L)<4095L)
  ;
  DWEL_projection_info=[DWEL_projection_info, $
    'star_r_Stats=('+strtrim(string(amin),2)+',' $
    +strtrim(string(amean),2)+',' $
    +strtrim(string(amax),2)+',' $
    +strtrim(string(asdev),2)+')' $
    ]
  ;
  image_statistics,star_r2,mask=mask,minimum=amin,maximum=amax,mean=amean,stddev=asdev
  scaler=4095.0/(amax-amin)
  star_r2=((round(scaler*(star_r2-amin))>0L)<4095L)
  ;
  DWEL_projection_info=[DWEL_projection_info, $
    'star_r2_Stats=('+strtrim(string(amin),2)+',' $
    +strtrim(string(amean),2)+',' $
    +strtrim(string(amax),2)+',' $
    +strtrim(string(asdev),2)+')' $
    ]
    
  ;set up the ENVI header for the output image and open in the
  ;available files list
    
  descrip=pname+' Projected BIL version of '+strtrim(DWEL_Cube_File,2)
  
  ;get output_file name without path
  last=strpos(DWEL_AT_File,path_sep(),/reverse_search)
  out_base=strtrim(strmid(DWEL_AT_File,last+1,strlen(DWEL_AT_File)-last-1),2)
  
  DWEL_projection_info=[DWEL_projection_info, $
    'Input_File='+strtrim(f_base,2),$
    'Output_projected_file='+strtrim(out_base,2)]
    
  bnames=band_name[band_pos]+'_'+ptype+'_Projected'
  
  envi_setup_head,fname=DWEL_AT_File,ns=ns_out,nl=nl_out,nb=nb_out,$
    xstart=0,ystart=0,$
    data_type=2, interleave=1, $
    descrip=descrip, wl=wl[band_pos], bnames=bnames,/write
    
  envi_open_file,DWEL_AT_File,r_fid=out_fid,/no_interactive_query,/no_realize
  
  ;write out the previous header records
  status=DWEL_put_headers(out_fid,DWEL_headers)
  
  envi_assign_header_value, fid=out_fid, keyword='DWEL_projection_info', $
    value=DWEL_projection_info
  envi_assign_header_value, fid=out_fid, keyword='DWEL_Andrieu_zenith', $
    value=DWEL_Andrieu_zenith
  envi_assign_header_value, fid=out_fid, keyword='DWEL_Andrieu_azimuth', $
    value=DWEL_Andrieu_azimuth
    
  envi_write_file_header, out_fid
  envi_file_mng,id=out_fid,/remove
  
  pos = 0b
  pos = where(mask eq 0, npos)
  if (npos gt 0) then begin
    maxwf[pos] = 0.0
    meanwf[pos] = 0.0
    accum[pos] = 0.0
    accum_r[pos] = 0.0
    accum_r2[pos] = 0.0
    star_r[pos] = 0.0
    star_r2[pos] = 0.0
  endif 

  ;now write out the extra information image
  n_base=strlen(DWEL_AT_File)
  n_dot=strpos(DWEL_AT_File,'.',/reverse_search)
  if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
    outextra[j]=DWEL_AT_File+'_extrainfo.img'
  endif else begin
    outextra=strmid(DWEL_AT_File,0,n_dot)+'_extrainfo.img'
  endelse
  ;Open output file
  text_err=0
  openw, ofile, outextra,/get_lun,error=text_err
  if (text_err ne 0) then begin
    print,'Error opening output file '+strtrim(outextra,2)
    err=12
    err_flag=1b
    goto, cleanup
  endif
  
  writeu,ofile,long(num_val)
  writeu,ofile,long(theta)
  writeu,ofile,long(phi)
  writeu,ofile,long(mask)
  writeu,ofile,long(maxwf)
  ;  writeu,ofile,fix(round(meanwf))
  writeu,ofile,long(accum)
  writeu,ofile,long(accum_r)
  writeu,ofile,long(accum_r2)
  writeu,ofile,long(star_r)
  writeu,ofile,long(star_r2)
  free_lun, ofile,/force
  
  descrip='Numbers and info for '+strtrim(DWEL_Cube_File,2)
  bnames=['Number Averaged','Zenith','Azimuth','Mask','Max','Mean','Mean_r','Mean_r2','Star_r','Star_RMS']
  ;  bnames=['Number Averaged','Zenith','Azimuth','Mask','Mean','Max']
  envi_setup_head,fname=outextra,ns=ns_out,nl=nl_out,nb=10,$
    xstart=0,ystart=0,$
    data_type=3, interleave=0, bnames=bnames, $
    descrip=descrip, /write
    
  envi_open_file,outextra,r_fid=anc_fid,/no_interactive_query,/no_realize
  
  ;write out the previous header records
  status=DWEL_put_headers(anc_fid,DWEL_headers)
  
  envi_assign_header_value, fid=anc_fid, keyword='DWEL_projection_info', $
    value=DWEL_projection_info
  envi_assign_header_value, fid=anc_fid, keyword='DWEL_Andrieu_zenith', $
    value=DWEL_Andrieu_zenith
  envi_assign_header_value, fid=anc_fid, keyword='DWEL_Andrieu_azimuth', $
    value=DWEL_Andrieu_azimuth
    
  envi_write_file_header, anc_fid
  envi_file_mng,id=anc_fid,/remove
  anc_fid=0b
  
  print,'Completed writing projected image - now for summary data'
  
  print,'Output File: '+strtrim(DWEL_AT_File,2)
  print,' '
  print,'DWEL_Projection_Info written to Output File Headers:'
  for index=0,n_elements(DWEL_projection_info)-1 do begin
    print,strtrim(DWEL_projection_info[index],2)
  endfor
  print,' '
  print,'*************************************'
  print,' '
  
  cleanup:
  
  num_val=0b
  theta=0b
  phi=0b
  mask=0b
  accum=0b
  accum_abs=0b
  pos=0b
  accum_r=0b
  accum_r2=0b
  star_r=0b
  star_r2=0b
  
  free_lun,inlun,/force
  free_lun,ofile,/force
  
  result=ptr_valid(p_stat)
  if (result) then begin
    ptr_free,p_stat
  endif
  
  result=ptr_valid(p_list[0])
  if (result) then begin
    ptr_free,p_list
  endif
  
  p_list=0b
  
  if (err_flag) then print,'Returning from dwel_cube2at with error'
  ;
  return
end
