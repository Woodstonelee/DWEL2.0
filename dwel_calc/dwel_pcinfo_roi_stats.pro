pro DWEL_pcinfo_roi_stats, pcinfo_file,roi_file,waveform_file,output_file,settings,err
  ;Get ROI stats from pcinfo file and waveform file
  ;Assumes hard targets such as panels and trunks
  compile_opt idl2
  
  tfile=98
  outfile=''
  ofile=99L
  fid=50
  wf_fid=55
  lib_fid=60
  err_flag=0b
  err=0
  
  ;clean up any fids which are no longer where they were!
  ;ENVI issue that is annoying and leads to confusion
  clean_envi_file_fids
  
  print,'entering pcinfo_roi_stats'
  
  ;do very protective checks any here returned with err=1
  pcinfo_file=strtrim(pcinfo_file,2)
  if (~file_test(pcinfo_file) or (strlen(pcinfo_file) le 0)) then begin
    print,'input pcinfo file does NOT exist!'
    print,'File entered='+strtrim(pcinfo_file,2)
    err_flag=1b
    err=1
    goto,cleanup
  endif
  ;
  roi_file=strtrim(roi_file,2)
  if (~file_test(roi_file) or (strlen(roi_file) le 0)) then begin
    print,'input roi file does NOT exist!'
    print,'File entered='+strtrim(roi_file,2)
    err_flag=1b
    err=1
    goto,cleanup
  endif
  ;
  waveform_file=strtrim(waveform_file,2)
  if (~file_test(waveform_file) or (strlen(waveform_file) le 0)) then begin
    print,'input waveform file does NOT exist!'
    print,'File entered='+strtrim(waveform_file,2)
    err_flag=1b
    err=1
    goto,cleanup
  endif
  ;
  save_waveforms=settings.save_waveforms
  
  ; Open pcinfo file
  envi_open_file, pcinfo_file, r_fid=fid,/no_realize,/no_interactive_query
  
  if (fid eq -1) then begin
    print,strtrim('Error opening input pcinfo file',2)
    print,'Input File: '+strtrim(pcinfo_file,2)
    err_flag=1b
    envi_file_mng,id=fid,/remove
    err=2
    goto, cleanup
  endif
  
  ; Get the file dimensions etc
  envi_file_query, fid, nb=nbands, nl=nlines, $
    ns=nsamples, bnames=bnames, data_type=dt, $
    interleave=ftype, dims=dims
    
  if (nbands ne 14) then begin
    print,'Input file not a pcinfo file - number of bands ne 14'
    err_flag=1b
    err=3
    goto,cleanup
  endif
  ;set the type of the file
  ft_nam='Unknown'
  case ftype of
    0: ft_nam='BSQ'
    1: ft_nam='BIL'
    2: ft_nam='BIP'
  endcase
  
  ;get path and name of pcinfo file as separate strings
  f_base=file_basename(pcinfo_file)
  f_path=file_dirname(pcinfo_file)
  
  print,'f_base='+strtrim(f_base,2)
  
  ;now get headers etc
  ;now get the DWEL headers that are present
  ;set up a base structure for the DWEL headers
  DWEL_headers={ $
    f_base:f_base $
    }
    
  ;find all of the DWEL headers in the pcinfo file as defined by FID
  status=DWEL_get_headers(fid,DWEL_headers)
  
  if (not status) then begin
    print,strtrim('Bad call to DWEL_get_headers! DWEL Header setup cancelled!',2)
    print,'Input File: '+strtrim(pcinfo_file,2)
    err_flag=1b
    err=4
    goto, cleanup
  endif
  
  if ((DWEL_headers.headers_present le 0s) or ~DWEL_headers.run_present) then begin
    print,strtrim('Input file is NOT a valid DWEL file!',2)
    print,'Input File: '+strtrim(pcinfo_file,2)
    envi_file_mng,id=fid,/remove
    err_flag=1b
    err=5
    goto, cleanup
  endif
  
  if (~DWEL_headers.base_present) then begin
    print,strtrim('Input file has NOT been basefixed!',2)
    print,'Input File: '+strtrim(pcinfo_file,2)
    envi_file_mng,id=fid,/remove
    err_flag=1b
    err=6
    goto, cleanup
  endif
  
  ;check for processing Level
  level1='File has been baseline fixed'
  level2='File has been saturation fixed'
  if (~(DWEL_headers.base_present) or ~(DWEL_headers.sat_present)) then begin
    if (~(DWEL_headers.base_present)) then level1='File has NOT been baseline fixed'
    if (~(DWEL_headers.sat_present)) then level2='File has NOT been saturation fixed'
    print,'File header indicates DWEL data has NOT been processed to the'
    print,'Level needed for Calibration. The Level needed is:'
    print,'DWEL data to be baseline fixed, saturation fixed and filtered.'
    print,'For the File '+f_base+':'
    print,level1
    print,level2
    envi_file_mng,id=fid,/remove
    err_flag=1b
    err=7
    goto, cleanup
  endif
  
  info=DWEL_headers.dwel_adaptation
  
  wavelength=1548
  
  ;now get the DWEL wavelength
  match = -1
  for i=0,n_elements(info)-1 do begin
    if (strmatch(info[i],'*Wavelength=*', /fold_case)) then match=i
  endfor
  if match ge 0 then begin
    text=strtrim(info[match],2)
    k=strpos(text,'=')
    wavelength=fix(strtrim(strmid(text,k+1,4),2))
  endif else begin
    print,strtrim('Input file has NOT got a wavelength record!',2)
    print,'Input File: '+strtrim(pcinfo_file,2)
    envi_file_mng,id=fid,/remove
    err_flag=1b
    err=18
    goto, cleanup
  endelse
  
  print,'pcinfo wavelength='+strtrim(string(wavelength),2)
  
  if (~DWEL_headers.ptcld_present) then begin
    print,strtrim('Input file has NOT been through point cloud processing!',2)
    print,'Input File: '+strtrim(pcinfo_file,2)
    envi_file_mng,id=fid,/remove
    err_flag=1b
    err=8
    goto, cleanup
  endif
  
  info=DWEL_headers.dwel_pointcloud_info
  
  ;now get the scaling information for the data
  match = -1
  for i=0,n_elements(info)-1 do begin
    if (strmatch(info[i],'*d_Stats=*', /fold_case)) then match=i
  endfor
  if match ge 0 then begin
    text=strtrim(info[match],2)
    k=strpos(text,'=')
    buf=strtrim(strmid(text,k+1),2)
    blen=strlen(buf)-2
    d_stats=float(strsplit(strmid(buf,1,blen),',',/extract))
  ;   print,d_stats
  endif else begin
    print,'d_Stats NOT present in pcinfo file!'
    print,'Input File: '+strtrim(pcinfo_file,2)
    envi_file_mng,id=fid,/remove
    err_flag=1b
    err=20
    goto, cleanup
  endelse
  match = -1
  for i=0,n_elements(info)-1 do begin
    if (strmatch(info[i],'*d0_Stats=*', /fold_case)) then match=i
  endfor
  if match ge 0 then begin
    text=strtrim(info[match],2)
    k=strpos(text,'=')
    buf=strtrim(strmid(text,k+1),2)
    blen=strlen(buf)-2
    d0_stats=float(strsplit(strmid(buf,1,blen),',',/extract))
  ;   print,d0_stats
  endif else begin
    print,'d0_Stats NOT present in pcinfo file!'
    print,'Input File: '+strtrim(pcinfo_file,2)
    envi_file_mng,id=fid,/remove
    err_flag=1b
    err=21
    goto, cleanup
  endelse
  match = -1
  for i=0,n_elements(info)-1 do begin
    if (strmatch(info[i],'*I_Stats=*', /fold_case)) then match=i
  endfor
  if match ge 0 then begin
    text=strtrim(info[match],2)
    k=strpos(text,'=')
    buf=strtrim(strmid(text,k+1),2)
    blen=strlen(buf)-2
    i_stats=float(strsplit(strmid(buf,1,blen),',',/extract))
  ;   print,i_stats
  endif else begin
    print,'I_Stats NOT present in pcinfo file!'
    print,'Input File: '+strtrim(pcinfo_file,2)
    envi_file_mng,id=fid,/remove
    err_flag=1b
    err=22
    goto, cleanup
  endelse
  
  d_gain=(d_stats[2]-d_stats[0])/4095.0d0
  d_off=d_stats[0]
  d0_gain=(d0_stats[2]-d0_stats[0])/4095.0d0
  d0_off=d0_stats[0]
  i_gain=(i_stats[2]-i_stats[0])/4095.0d0
  i_off=i_stats[0]
  
  ;now the waveform file
  envi_open_file, waveform_file, r_fid=wf_fid,/no_realize,/no_interactive_query
  
  if (wf_fid eq -1) then begin
    print,strtrim('Error opening input waveform file',2)
    print,'Input File: '+strtrim(waveform_file,2)
    err_flag=1b
    envi_file_mng,id=fid,/remove
    err=9
    goto, cleanup
  endif
  
  ; Get the file dimensions etc
  envi_file_query, wf_fid, nb=nb, nl=nl, $
    ns=ns, data_type=type,wl=range
    
  if ((nl ne nlines) or (ns ne nsamples)) then begin
    print,'Input waveform file does NOT conform with pcinfo'
    print,'Input PcInfo File='+strtrim(pcinfo_file,2)
    print,'Input Waveform File='+strtrim(waveform_file,2)
    envi_file_mng,id=fid,/remove
    envi_file_mng,id=wf_fid,/remove
    err_flag=1b
    err=10
    goto,cleanup
  endif
  
  ;get number of bytes in input waveform file
  nbytes=dt2nb(type)
  
  w_base=file_basename(waveform_file)
  print,'w_base='+strtrim(w_base,2)
  
  ;now get headers etc
  ;now get the DWEL headers that are present
  ;set up a base structure for the DWEL headers
  wf_headers={ $
    f_base:w_base $
    }
    
  ;find all of the DWEL headers in the pcinfo file as defined by FID
  status=DWEL_get_headers(wf_fid,wf_headers)
  
  if (not status) then begin
    print,strtrim('Bad call to DWEL_get_headers! DWEL Header setup cancelled!',2)
    print,'Input File: '+strtrim(waveform_file,2)
    envi_file_mng,id=fid,/remove
    envi_file_mng,id=wf_fid,/remove
    err_flag=1b
    err=11
    goto, cleanup
  endif
  
  if ((wf_headers.headers_present le 0s) or ~wf_headers.run_present) then begin
    print,strtrim('Input waveform file is NOT a valid DWEL file!',2)
    print,'Input File: '+strtrim(waveform_file,2)
    envi_file_mng,id=fid,/remove
    envi_file_mng,id=wf_fid,/remove
    err_flag=1b
    err=12
    goto, cleanup
  endif
  
  if (~wf_headers.base_present) then begin
    print,strtrim('Input waveform file has NOT been basefixed!',2)
    print,'Input File: '+strtrim(waveform_file,2)
    envi_file_mng,id=fid,/remove
    envi_file_mng,id=wf_fid,/remove
    err_flag=1b
    err=13
    goto, cleanup
  endif
  
  ;check for waveform file processing Level
  level1='File has been baseline fixed'
  level2='File has been saturation fixed'
  if (~(wf_headers.base_present) or ~(wf_headers.sat_present)) then begin
    if (~(wf_headers.base_present)) then level1='File has NOT been baseline fixed'
    if (~(wf_headers.sat_present)) then level2='File has NOT been saturation fixed'
    print,'File header indicates waveform data has NOT been processed to the'
    print,'Level needed for Calibration. The Level needed is:'
    print,'DWEL data to be baseline fixed, saturation fixed and filtered.'
    print,'For the File '+w_base+':'
    print,level1
    print,level2
    envi_file_mng,id=fid,/remove
    envi_file_mng,id=wf_fid,/remove
    err_flag=1b
    err=14
    goto, cleanup
  endif
  
  info=wf_headers.dwel_adaptation
  
  wf_wavelength=1548
  ;now get the DWEL wavelength
  match = -1
  for i=0,n_elements(info)-1 do begin
    if (strmatch(info[i],'*Wavelength=*', /fold_case)) then match=i
  endfor
  if match ge 0 then begin
    text=strtrim(info[match],2)
    k=strpos(text,'=')
    wf_wavelength=fix(strtrim(strmid(text,k+1,4),2))
  endif else begin
    print,strtrim('Input waveform file has NOT got a wavelength record!',2)
    print,'Input File: '+strtrim(waveform_file,2)
    envi_file_mng,id=fid,/remove
    envi_file_mng,id=wf_fid,/remove
    err_flag=1b
    err=19
    goto, cleanup
  endelse
  
  print,'waveform wavelength='+strtrim(string(wf_wavelength),2)
  
  if (wf_wavelength ne wavelength) then begin
    print,strtrim('Input waveform Wavelength NOT the same as the PcInfo!',2)
    print,'Input PcInfo File='+strtrim(pcinfo_file,2)
    print,'Input Waveform File='+strtrim(waveform_file,2)
    envi_file_mng,id=fid,/remove
    envi_file_mng,id=wf_fid,/remove
    err_flag=1b
    err=19
    goto, cleanup
  endif
  
  ;set up and open the output info file
  ;see if the output file exists & remove if it does!
  if(file_test(output_file)) then begin
    fids=envi_get_file_ids()
    if(fids[0] eq -1) then begin
      file_delete, output_file,/quiet
      print,'old output file deleted'
    endif else begin
      for i=0,n_elements(fids)-1 do begin
        envi_file_query,fids[i],fname=tname
        if (strtrim(strlowcase(output_file),2) eq $
          strtrim(strlowcase(tname),2)) then begin
          envi_file_mng,id=fids[i],/remove
          print,'old output file removed from ENVI'
        endif
      endfor
      file_delete, output_file,/quiet
      print,'old output file deleted'
    endelse
  endif
  text_err=0
  openw, tfile, output_file,/get_lun,error=text_err
  if (text_err ne 0) then begin
    print,'Error opening output file in point cloud!!'
    print,'File Name =',strtrim(output_file,2)
    print,'text_err=',text_err
    print,'Error Type =',strtrim(string(!ERROR_STATE.MSG),2)
    print,'Sys_Error Type =',strtrim(string(!ERROR_STATE.SYS_MSG),2)
    err_flag=1b
    err=15
    goto, cleanup
  endif
  
  time_date=strtrim(systime(),2)
  
  printf,tfile,strtrim('[DWEL PcInfo ROI Stats Data]',2)
  printf,tfile,strtrim('Run made at: '+time_date,2)
  printf,tfile,strtrim('Input PcInfo File='+pcinfo_file,2)
  printf,tfile,strtrim('Input ROI File='+roi_file,2)
  printf,tfile,strtrim('Input Waveform File='+waveform_file,2)
  flush,tfile
  
  ;now set up the rois
  envi_restore_rois,roi_file
  cal_roi_ids=envi_get_roi_ids(fid=roi_fid,roi_colors=roi_colors,roi_names=roi_names)
  if (cal_roi_ids[0] eq -1) then begin
    print,'cal_roi_ids=',cal_roi_ids
    print,'rois not set up, exit!'
    err_flag=1b
    err=16
    goto, cleanup
  endif
  help,roi_colors
  num_rois=n_elements(cal_roi_ids)
  print,'number of ROIs='+strtrim(string(num_rois),2)
  print,'roi_id,roi_color,roi_name'
  for i=0,num_rois-1 do begin
    temp=reform(roi_colors[*,i])
    buf=string(temp,format='(i3,",",i3,",",i3)')
    buf='['+strtrim(buf,2)+']'
    print,strtrim(string(cal_roi_ids[i]),2),',',strtrim(buf,2),',',strtrim(roi_names[i],2)
  endfor
  
  printf,tfile,'number of ROIs='+strtrim(string(num_rois),2)
  printf,tfile,'roi_id,roi_color,roi_name'
  for i=0,num_rois-1 do begin
    temp=reform(roi_colors[*,i])
    buf=string(temp,format='(i3,",",i3,",",i3)')
    buf='['+strtrim(buf,2)+']'
    printf,tfile,strtrim(string(cal_roi_ids[i]),2),',',strtrim(buf,2),',',strtrim(roi_names[i],2)
  endfor
  flush,tfile
  
  ;get total_d,total_d0,total_I,mean_range,zenith,azimuth,residual,z
  anc_pos=indgen(nbands)
  
  ;now do the business
  printf,tfile,'ROI information'
  printf,tfile,'ROI_Name,nsamp,Mean_d,Mean_d0,mean_I,mean_Range,Mean_zenith,Mean_Azimuth,Mean_residual,Mean_z'
  print,'ROI_Name,nsamp,Mean_d,Mean_d0,mean_I,mean_Range,Mean_zenith,Mean_Azimuth,Mean_residual,Mean_z'
  for i=0,num_rois-1 do begin
    result=float(envi_get_roi_data(cal_roi_ids[i],fid=fid,pos=anc_pos))
    nsamp=(size(result))[2]
    pos_valid=where((result[0,*] eq 1) and (result[11,*] ne 1) and (result[12,*] eq 1),n_valid)
    if (n_valid gt 0) then begin
      nsamp=n_valid
      result=result[*,pos_valid]
      mean_roi=total(result,2)/float(nsamp)
      roi_sq=total(result^2,2)/float(nsamp)
      min_roi=min(result,dimension=2,max=max_roi)
      ;
      outstring=strtrim(roi_names[i],2)+','+$
        strtrim(string(nsamp,format='(i10)'),2)+','+$
        strtrim(string(d_gain*mean_roi[1]+d_off,format='(f10.2)'),2)+','+$
        strtrim(string(d0_gain*mean_roi[2]+d0_off,format='(f10.2)'),2)+','+$
        strtrim(string(i_gain*mean_roi[3]+i_off,format='(f10.2)'),2)+','+$
        strtrim(string(mean_roi[6]/100.0,format='(f10.2)'),2)+','+$
        strtrim(string(mean_roi[8]/100.0,format='(f10.2)'),2)+','+$
        strtrim(string(mean_roi[9]/100.0,format='(f10.2)'),2)+','+$
        strtrim(string(mean_roi[10],format='(f10.2)'),2)+','+$
        strtrim(string(mean_roi[13]/100.0,format='(f10.2)'),2)
      print,strtrim(outstring,2)
      printf,tfile,strtrim(outstring,2)
    endif else begin
      outstring=strtrim(roi_names[i],2)+',0,0.0,0.0,0.0,0,0.0,0.0,0.0,0.0'
      print,strtrim(outstring,2)
      printf,tfile,strtrim(outstring,2)
    endelse
    result=0b
  endfor
  flush,tfile
  result=0b
  
  ;================================================================
  ;Now if desired open up the spectral library file
  ;for the means after removing bad FWHM data
  if (save_waveforms) then begin
    n_base=strlen(waveform_file)
    n_dot=strpos(waveform_file,'.',/reverse_search)
    if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
      outlib=waveform_file+'_waveform_roi_means.sli'
    endif else begin
      outlib=strmid(waveform_file,0,n_dot)+'_waveform_roi_means.sli'
    endelse
    ;see if the output file exists and delete if it does
    if(file_test(outlib)) then begin
      fids=envi_get_file_ids()
      if(fids[0] eq -1) then begin
        file_delete, outlib
      endif else begin
        for i=0,n_elements(fids)-1 do begin
          envi_file_query,fids[i],fname=tname
          if (strtrim(strlowcase(outlib),2) eq $
            strtrim(strlowcase(tname),2)) then begin
            envi_file_mng,id=fids[i],/remove
          endif
        endfor
        file_delete, outlib
      endelse
    endif
    ;Open file for Spectral Library
    text_err=0
    openw, ofile, outlib,/get_lun,error=text_err
    if (text_err ne 0) then begin
      print,'Error opening output spectral library file '+strtrim(outfile,2)
      print,'Error Name is: '+strtrim(!error_state.name,2)
      goto, cleanup
    endif
  endif
  ;================================================================
  
  range_pos=where(range gt 0.0,n_rpos)
  wf_pos=indgen(nb)
  w_step=range[1]-range[0]
  
  ;now process the waveform information
  printf,tfile,'ROI Waveform information'
  printf,tfile,'ROI_Name,nsamp,Max_d,Range,FWHM'
  print,'ROI_Name,nsamp,Max_d,Range,FWHM'
  for i=0,num_rois-1 do begin
    result=float(envi_get_roi_data(cal_roi_ids[i],fid=fid,pos=anc_pos))
    nsamp=(size(result))[2]
    pos_valid=where((result[0,*] eq 1) and (result[11,*] ne 1) and (result[12,*] eq 1),n_valid)
    if (n_valid gt 0) then begin
      result=0b
      nsamp=n_valid
      result=float(envi_get_roi_data(cal_roi_ids[i],fid=wf_fid,pos=wf_pos))
      result=result[*,pos_valid]
      mean_roi=total(result,2)/float(nsamp)
      ;    roi_sq=total(result^2,2)/float(nsamp)
      if (save_waveforms) then writeu,ofile,mean_roi
      max_val=max(float(mean_roi[range_pos]),ipos)
      fwhm=w_step*total(float(mean_roi[range_pos]))/max_val
      ;
      outstring=strtrim(roi_names[i],2)+','+$
        strtrim(string(nsamp,format='(i10)'),2)+','+$
        strtrim(string(max_val,format='(f10.2)'),2)+','+$
        strtrim(string(range[range_pos[ipos]],format='(f10.2)'),2)+','+$
        strtrim(string(fwhm,format='(f10.2)'),2)
      print,strtrim(outstring,2)
      printf,tfile,strtrim(outstring,2)
    endif else begin
      result=0b
      mean_roi=fltarr(nb)
      if (save_waveforms) then writeu,ofile,mean_roi
      outstring=strtrim(roi_names[i],2)+',0,0.0,0.0,0.0'
      print,strtrim(outstring,2)
      printf,tfile,strtrim(outstring,2)
    endelse
    result=0b
    pos_valid=0b
    outstring=''
    nsamp=0b
  endfor
  flush,tfile
  result=0b
  
  ;now if rois being saved set up header etc
  if (save_waveforms) then begin
    free_lun,ofile,/force
    descrip='Mean ROI Spectra for '+strtrim(waveform_file,2)
    bnames=['Mean ROI Spectra']
    envi_setup_head,fname=outlib,$
      ns=nb, nl=num_rois,nb=1,file_type=4,interleave=0, $
      data_type=4,wl=range,bnames=bnames, $
      spec_names=roi_names,descrip=descrip, $
      zplot_titles=['Range(m)','Value'],$
      /write
    envi_open_file,outlib,r_fid=lib_fid,/no_realize,/no_interactive_query
    ;write out the previous header records
    status=dwel_put_headers(lib_fid,DWEL_headers)
    envi_file_mng,id=lib_fid,/remove
  endif
  
  envi_file_mng,id=fid,/remove
  envi_file_mng,id=wf_fid,/remove
  
  free_lun,tfile,/force
  
  cleanup:
  help,/memory
  ;clean up any rois
  roi_ids = envi_get_roi_ids()
  if (roi_ids[0] ge 0)then begin
    print,'number of rois removed=',n_elements(roi_ids)
    envi_delete_rois,/all
  endif
  ;free units and pointers
  envi_file_mng,id=fid,/remove
  envi_file_mng,id=wf_fid,/remove
  envi_file_mng,id=lib_fid,/remove
  free_lun,tfile,/force
  free_lun,ofile,/force
  
  if (err_flag) then begin
    print,'pcinfo_roi_stats finished with error'
  endif else begin
    print,'pcinfo_roi_stats finished without error'
  endelse
  heap_gc,/verbose
  return
;
end
