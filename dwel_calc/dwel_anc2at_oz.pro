; Simple AT projection of unprocessed data cube.
; Arrange pixels in array of zenith and azimuth angle instead of shot number and scan number

pro DWEL_val_block
  pi_val=4.0*atan(1.0)
  rad=pi_val/180.0
  deg=1.0/rad
  eta=1.0e-7
end

; Max_Zenith_Angle: in unit of degree
; output_resolution: in unit of mrad
pro DWEL_anc2at_oz, DWEL_Anc_File, DWEL_AT_File, Max_Zenith_Angle, output_resolution, zen_tweak, err

  compile_opt idl2
;  envi, /restore_base_save_files
;  envi_batch_init, /no_status_window

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
  err=0
  ;;;;;;;;;;;;;;;;;;;;;;;;
  
  envi_open_file, DWEL_Anc_File,r_fid=anc_fid,/no_interactive_query,/no_realize
  if (anc_fid eq -1) then begin
    print,'Processing stopped! Error opening ancillary data file '+strtrim(DWEL_Anc_File,2)
    err=1
    goto, cleanup
  endif
  
  ;get the input image dimensions and other info
  envi_file_query, anc_fid, ns=Nshots, nl=Nscans, nb=nb_anc, $
                   data_type=type, file_type=ftype, dims=dims

  samples=Nshots
  lines=Nscans
  bands=nb_anc
  
  print,'ancfile type=',ftype
  
  band_pos=indgen(bands)
  
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
  last=strpos(DWEL_Anc_File,path_sep(),/reverse_search)
  f_path=strmid(DWEL_Anc_File,0,last+1)
  f_base=strtrim(strmid(DWEL_Anc_File,last+1,strlen(DWEL_Anc_File)-last-1),2)
  ;now get the EVI headers that are present
  ;set up a base structure for the EVI headers
  DWEL_headers={ $
       f_base:f_base $
       }
  
  ;find all of the EVI headers in the hdr file as defined by FID
  status=DWEL_get_headers(anc_fid,DWEL_headers)
  
  if (not status) then begin
    print, 'Processing stopped! Bad FID in EVI get_headers on input file!'
    err=2
    goto, cleanup
  endif
  
  if (DWEL_headers.headers_present le 0s or not DWEL_headers.run_present) then begin
    print,'Processing stopped! Input File is NOT a valid EVI Cube file!'
    err=3
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

  print,'Input information from input files complete'

  ;get the mask
  Mask_all=bytarr(Nshots,Nscans)+1b
  dims=[-1,0,Nshots-1,0,Nscans-1]
  Mask_all=byte(envi_get_data(fid=anc_fid,dims=dims,pos=6))
  ;; get the waveform maximum image
  wfmax = long(envi_get_data(fid=anc_fid,dims=dims,pos=5))
  
  ;get encoder positions!
  ShotZen=fltarr(Nshots,Nscans)
  ShotAZim=fltarr(Nshots,Nscans)
  ShotZen=float(envi_get_data(fid=anc_fid,dims=dims,pos=2))
  ShotAZim=float(envi_get_data(fid=anc_fid,dims=dims,pos=3))
  
  envi_file_mng,id=anc_fid,/remove

sel=0
bad=50
bad_end=10
RotEnc_max=524288.0d0
Rot_med=double(median(ShotAzim,dimension=1))
;
Rot_med=reform(Rot_med[bad:Nscans-1])
pos=where(Rot_med gt 0,npos)
if (npos gt 0) then Rot_med=reform(Rot_med[pos])
numval=n_elements(Rot_med)

Rot_med_compl=Rot_med-RotEnc_Max/2.0d0
posc=where(Rot_med_compl lt 0.0,nposc)
if (nposc gt 0) then Rot_med_compl[posc]=Rot_med_compl[posc]+RotEnc_Max
posc=0b

valend=Rot_med[numval-bad_end-1]
posc=where(Rot_med_compl lt valend,nposc)
if (nposc gt 0) then begin
  sel=max([bad,pos[posc[0]]+bad])
endif else begin
  sel=bad
endelse
posc=0b
pos=0b
print,'final sel=',sel
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

  status = DWEL_set_theta_phi_oz(p_stat,zen_tweak)
  
  ;put the results into the local arrays
  ShotZen=(*p_stat).ShotZen
  ShotAzim=(*p_stat).ShotAzim
  ptr_free, p_stat
  
  ;set the mask for angles outside range
  pos=where(Shotzen gt Max_Zenith_Angle,npos)
  print,'npos angles above max=',npos
  
  if (npos gt 0) then Mask_all[pos]=0
  pos=0b
  
  ;now introduce the new mask to ensure minimum wrap
  mask_all[*,0:sel]=0
  mask_all[*,(nscans-bad_end):(nscans-1)]=0
  
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
  print,'sampling rate=',srate

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
  
  ;print, 'Andrieu_transpose'
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
  
  ;now loop over the cells and set up pointers to the original image cells
  ;involved in the output spatial cells
  for i=0, nl_out-1 do begin
    y=float(i)*h
    y_min=y-h2
    y_max=y+h2
    pos_y=where(((y_proj ge y_min) and (y_proj le y_max) and Mask_all),count_y)

    if (count_y gt 0) then begin
      for j=0,ns_out-1 do begin
        x=float(j)*h
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
            pos_ind=[[pos_ind],[temp]] ; pos_ind is n*3 array, n is the number of shots in the line i of the projection, first and second columns are the pixel location in the original data cube, the third line is the column j of the projection. 
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
        theta[j,i]=round(100.0*th)
        phi[j,i]=round(100.0*ph)
      endfor
      if (total(reform(num_val[*,i])) gt 0) then begin
        temp=sort(reform(pos_ind[1,*])) ; sort the pos_ind by column, the shot number/zenith
        pos_ind=pos_ind[*,temp]
      endif
    endif else begin
      for j=0,ns_out-1 do begin
        x=float(j)*h
        num_val[j,i]=0L
        mask[j,i]=0L
        status=t_Andrieu_xy2tp(1,x,y,th,ph)
        A_theta[j]=th
        A_phi[i]=ph
        theta[j,i]=round(100.0*th)
        phi[j,i]=round(100.0*ph)
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

  DWEL_Andrieu_zenith=a_zen
  DWEL_Andrieu_azimuth=a_azm

  a_zen=0b
  a_azm=0b
  A_theta=0b
  A_phi=0b
  
  if (type le 1) then begin
    scaler=100.0d0
  endif else begin
    scaler=1.0d0
  endelse
  
  DWEL_Anc2AT_info=[ $
  'Program='+'DWEL Ancillary to AT projection Qlook',$
  'Processing Date Time='+strtrim(systime(),2),$
  'Projection type='+ptype,$
  'Projection name='+pname,$
  'Beam_Divergence_(mrad)='+strtrim(string(beam_div,format='(f14.3)'),2),$
  'Scan_Step_(mrad)='+strtrim(string(scan_step,format='(f14.3)'),2),$
  'Sampling_Ratio='+strtrim(string(sampling_ratio,format='(f14.3)'),2),$
  'output_resolution_(mrad)='+strtrim(string(ifov_x,format='(f10.2)'),2),$
  'max_zenith_angle_(deg)='+strtrim(string(!radeg*t_max,format='(f10.2)'),2),$
  'Zen_tweak_(enc)='+strtrim(string(zen_tweak),2),$
  'Mean image scale='+strtrim(string(scale,format='(f10.2)'),2),$
  'Output scale='+strtrim(string(scaler,format='(f10.2)'),2) $
  ]
  DWEL_Anc2AT_info=strtrim(DWEL_Anc2AT_info,2)
  
  ;all ready to go ... so get output file name[s]
  output_envi:
    
  print,'pre-processing done - projecting the image!'
    
  ;set up the arrays for data and output
  data=make_array(Nshots,1,/double)
  temp=make_array(Nshots,1,/double)
  maxwf = make_array(ns_out, nl_out, /double)

  num_avg=make_array(ns_out,nl_out,/long)
  
  ;do the processing over the output tiles
  
  for k=0L, nl_out-1 do begin
    current=-1
    ;count=0L
    temp=make_array(ns_out,1,/double)
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
          data = double(wfmax[*, lin])
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
    
    maxwf[*,k] = temp
    temp=0b
  endfor
  
  data=0b
  temp=0b
  
  if (total(abs(num_avg-num_val)) gt 0) then begin
    print,'numbers of averaged cells do NOT agree!'
    print,'total counted in first loop=',total(num_val)
    print,'total counted in second loop=',total(num_avg)
  endif
  
  image_statistics,maxwf,mask=mask,minimum=amin,maximum=amax,mean=amean,stddev=asdev
  scaler=4095.0/(amax-amin)
  maxwf_out=scaler*(maxwf-amin)
  
  pos=0b
  pos=where(mask eq 0,npos)
  if (npos gt 0) then begin
    maxwf_out[pos]=0.0
  endif
  
  print,'amin,amax,amean=',amin,amax,amean
  
  pos=0b
    
  ;now write out the extra information image
  ;; set up the file name of the extra information image
  n_base=strlen(DWEL_AT_File)
  n_dot=strpos(DWEL_AT_File,'.',/reverse_search)
;  if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
;    outextra[j]=DWEL_AT_File+'_extrainfo.img'
;  endif else begin
;    outextra=strmid(DWEL_AT_File,0,n_dot)+'_extrainfo.img'
;  endelse
  outextra=DWEL_AT_File
  ;Open output file
  text_err=0
  openw, ofile, outextra,/get_lun,error=text_err
  if (text_err ne 0) then begin
    print,'Error opening output file '+strtrim(outextra,2)
    err=4
    goto, cleanup
  endif
  
  writeu,ofile,long(num_val)
  writeu,ofile,long(theta)
  writeu,ofile,long(phi)
  writeu,ofile,long(mask)
;;  writeu,ofile,fix(round(meanwf))
;  writeu,ofile,round(maxwf_out)
  writeu,ofile,round(maxwf)
  free_lun, ofile,/force
  
  descrip='Numbers and info for '+strtrim(DWEL_Anc_File,2)
;  bnames=['Number Averaged','Zenith','Azimuth','Mask','Mean','Mean_r','Mean_r2','Star_r','Star_RMS']
  bnames=['Number Averaged','Zenith','Azimuth','Mask','Max']
  envi_setup_head,fname=outextra,ns=ns_out,nl=nl_out,nb=5,$
  xstart=0,ystart=0,$
  data_type=3, interleave=0, bnames=bnames, $
  descrip=descrip, /write
  
  envi_open_file,outextra,r_fid=anc_fid,/no_interactive_query,/no_realize
  
  ;write out the previous header records
  status=DWEL_put_headers(anc_fid,DWEL_headers)
  
  envi_assign_header_value, fid=anc_fid, keyword='DWEL_Anc2AT_info', $
      value=DWEL_Anc2AT_info
  
  envi_write_file_header, anc_fid
  envi_file_mng,id=anc_fid,/remove
  anc_fid=0b

  print,'Completed writing projected image - now for summary data'
  
  print,'Output File: '+strtrim(DWEL_AT_File,2)
  print,' '
  print,'DWEL_Anc2AT_info written to Output File Headers:'
  for index=0,n_elements(DWEL_Anc2AT_info)-1 do begin
    print,strtrim(DWEL_Anc2AT_info[index],2)
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
  
  result=ptr_valid(p_stat)
  if (result) then begin
    ptr_free,p_stat
  endif
  
  result=ptr_valid(p_list[0])
  if (result) then begin
    ptr_free,p_list
  endif
  p_list=0b
  
  if (err gt 0) then print,'Error called from dwel_anc2at'
  
end