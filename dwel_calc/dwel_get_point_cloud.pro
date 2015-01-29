;======================================================================
pro dwel_apply_ptcl_filter, p, pb_stats, pb_meta, pb_info, error=error 

  compile_opt idl2

  ;set data
  zero=0.0
  error=0
  point_file=''
  mdata_file=''
  mfile=99
  tfile=101
  tempfile=103
  bfile=105
  rfile=107
  pfile=108
  ppfile=109
  fid=30
  ovdebug=0b
  inlun=105
  
  ;set local parameters
  num_test=2
  test_num=20
  range_extent=5.0 ; unit, meter, half extent for pulse integral calculation
  noise_hits=10
  
  ;set height and width sieve
  zhigh=(*pb_meta).zhigh
  zlow=(*pb_meta).zlow
  xmin=(*pb_meta).xmin
  xmax=(*pb_meta).xmax
  ymin=(*pb_meta).ymin
  ymax=(*pb_meta).ymax
  save_br=(*pb_meta).save_br
  save_pfilt=(*pb_meta).save_pfilt
  DWEL_AppRefl=(*pb_meta).DWEL_AppRefl
  
  if (save_br) then ovdebug=1b else ovdebug=0b
  
  ;set other stuff from the meta_data
  DWEL_num=long((*pb_meta).run_number)
  threshold=(*pb_meta).threshold
  b_thresh=(*pb_meta).b_thresh
  r_thresh=(*pb_meta).r_thresh ; 0.175
  sieve_thresh=(*pb_meta).sieve_thresh
  wavelength=(*pb_meta).wavelength
  fneg=(*pb_meta).fneg
  h1=(*pb_meta).h1
  h2=(*pb_meta).h2
  i_scale=(*pb_meta).i_scale
  infile=(*pb_stats).infile
  if ((*pb_meta).Zero_Hit_Option) then add='_addgap' else add=''
  if ((*pb_meta).Add_DWEL) then add=add+'_adddwel' else add=add+''

  laser_man = (*pb_meta).laser_man

  cvthresh=threshold/512.0
  
  ; Open input and ancillary files
  envi_open_file,infile,r_fid=fid,/no_realize,/no_interactive_query
  if (fid eq -1) then begin
    print,strtrim('Error opening input DWEL file',2)
    print,'Input File: '+strtrim(infile,2)
    error=1
    goto,cleanup
  endif
  
  ; Get the file dimensions etc
  envi_file_query, fid,fname=infile,nb=nbands,nl=nlines, $
    ns=nsamples,bnames=bnames,descrip=descrip,wl=wl,data_type=dt
    
  ;get number of bytes in input file data
  nbytes=dt2nb(dt)

  if (save_br or save_pfilt) then begin
    f_base=file_basename(infile)
  ;set up a base structure for the DWEL headers
    DWEL_headers={ $
         f_base:f_base $
       }
  ;find all of the DWEL headers in the hdr file as defined by FID
    status=dwel_get_headers(fid,DWEL_headers)
  endif

  ;now remove the fid
  envi_file_mng,id=fid,/remove
  
  ;now set up the point and metadata files
  point_file=strtrim((*pb_stats).outfile,2)
  n_base=strlen(point_file)
  n_dot=strpos(point_file,'.',/reverse_search)
  if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
    point_file=strtrim(point_file,2)+'_points.txt'
    mdata_file=strtrim(point_file,2)+'_metadata.txt'
  endif else begin
    point_file=strtrim(strmid(point_file,0,n_dot),2)+'_points.txt'
    mdata_file=strtrim(strmid(point_file,0,n_dot),2)+'_metadata.txt'
  endelse
  
  ;see if the points file exists & remove if it does!
  if(file_test(point_file)) then begin
    fids=envi_get_file_ids()
    if(fids[0] eq -1) then begin
      file_delete, point_file,/quiet
      print,'old points file deleted'
    endif else begin
      for i=0,n_elements(fids)-1 do begin
        envi_file_query,fids[i],fname=tname
        if (strtrim(strlowcase(point_file),2) eq $
          strtrim(strlowcase(tname),2)) then begin
          envi_file_mng,id=fids[i],/remove
          print,'old points file removed from ENVI'
        endif
      endfor
      file_delete, point_file,/quiet
      print,'old points file deleted'
    endelse
  endif
  text_err=0
  openw, tfile, point_file,/get_lun,error=text_err
  if (text_err ne 0) then begin
    print,'Error opening points file in point cloud!!'
    print,'File Name =',strtrim(point_file,2)
    print,'text_err=',text_err
    print,'Error Type =',strtrim(string(!ERROR_STATE.MSG),2)
    print,'Sys_Error Type =',strtrim(string(!ERROR_STATE.SYS_MSG),2)
    error=2
    goto, cleanup
  endif
  
  time_date=strtrim(systime(),2)
  
  printf,tfile,strtrim('[DWEL Point Cloud Data]',2)
  printf,tfile,strtrim('Run made at: '+time_date,2)
  printf,tfile,strtrim('X,Y,Z,d_I,Return_Number,Number_of_Returns,Shot_Number,Run_Number,range,theta,phi,rk,Sample,Line,Band,FWHM',2)
  flush,tfile
  
  ;see if the metadata file exists & remove if it does!
  if(file_test(mdata_file)) then begin
    fids=envi_get_file_ids()
    if(fids[0] eq -1) then begin
      file_delete, mdata_file,/quiet
      print,'old meta data file deleted'
    endif else begin
      for i=0,n_elements(fids)-1 do begin
        envi_file_query,fids[i],fname=tname
        if (strtrim(strlowcase(mdata_file),2) eq $
          strtrim(strlowcase(tname),2)) then begin
          envi_file_mng,id=fids[i],/remove
          print,'old meta data file removed from ENVI'
        endif
      endfor
      file_delete, mdata_file,/quiet
      print,'old meta data file deleted'
    endelse
  endif
  text_err=0
  openw, mfile, mdata_file,/get_lun,error=text_err
  if (text_err ne 0) then begin
    print,'Error opening metadata file in point cloud!!'
    print,'File Name =',strtrim(mdata_file,2)
    print,'text_err=',text_err
    print,'Error Type =',strtrim(string(!ERROR_STATE.MSG),2)
    print,'Sys_Error Type =',strtrim(string(!ERROR_STATE.SYS_MSG),2)
    error=3
    goto, cleanup
  endif
  
  (*pb_meta).Processing_Date_Time=time_date
  
  printf,mfile,strtrim('[DWEL Point Cloud Metadata]',2)
  
  ;write out the metadata that you have set
  printf,mfile,'Processing_Date_Time='+strtrim((*pb_meta).Processing_Date_Time,2)
  printf,mfile,'Run_Number='+strtrim(string((*pb_meta).Run_Number),2)
  printf,mfile,'Description='+strtrim((*pb_meta).Description,2)
  printf,mfile,'Lasers='+strtrim((*pb_meta).laser_man,2)
  printf,mfile,'Wire_Flag='+strtrim(string((*pb_meta).wire_flag,format='(f10)'),2)
  printf,mfile,'Input_Path='+strtrim((*pb_meta).Input_Path,2)
  printf,mfile,'Input_File='+strtrim((*pb_meta).Input_File,2)
  printf,mfile,'Acquisition_Date_Time=',strtrim((*pb_meta).Acquisition_Date_Time,2)
  printf,mfile,'Ancillary_Path='+strtrim((*pb_meta).Ancillary_Path,2)
  printf,mfile,'Ancillary_File='+strtrim((*pb_meta).Ancillary_File,2)
  printf,mfile,'Projection='+strtrim((*pb_meta).Projection,2)
  printf,mfile,'DWEL_AppRefl='+strtrim(string((*pb_meta).DWEL_AppRefl),2)
  printf,mfile,'DWEL_Height(m)='+strtrim(string((*pb_meta).DWEL_Height,format='(f10.2)'),2)
  printf,mfile,'DWEL_Az_North='+strtrim(string((*pb_meta).DWEL_Az_North,format='(f10.2)'),2)
  printf,mfile,'Max_Zenith_Angle='+strtrim(string((*pb_meta).Max_Zenith_Angle,format='(f10.2)'),2)
  printf,mfile,'Range_Step(m)='+strtrim(string((*pb_meta).Range_Step,format='(f10.4)'),2)
  printf,mfile,'Threshold='+strtrim(string((*pb_meta).Threshold,format='(f10.4)'),2)
  printf,mfile,'b_thresh='+strtrim(string((*pb_meta).b_thresh,format='(f10.4)'),2) 
  printf,mfile,'sieve_thresh='+strtrim(string((*pb_meta).sieve_thresh,format='(f10.4)'),2)
  printf,mfile,'r_thresh='+strtrim(string((*pb_meta).r_thresh,format='(f10.4)'),2)
  printf,mfile,'Fneg='+strtrim(string((*pb_meta).Fneg,format='(f10.4)'),2)
  printf,mfile,'h1='+strtrim(string((*pb_meta).h1,format='(i10)'),2)
  printf,mfile,'h2='+strtrim(string((*pb_meta).h2,format='(i10)'),2)
  printf,mfile,'Wavelength='+strtrim(string((*pb_meta).wavelength,format='(i10)'),2)
  printf,mfile,'Zero_Hit_Option='+strtrim(string((*pb_meta).Zero_Hit_Option,format='(i10)'),2)
  printf,mfile,'Add_DWEL='+strtrim(string((*pb_meta).Add_DWEL,format='(i10)'),2)
  printf,mfile,'X_scale='+strtrim(string((*pb_meta).X_scale,format='(f10.3)'),2)
  printf,mfile,'Y_scale='+strtrim(string((*pb_meta).Y_scale,format='(f10.3)'),2)
  printf,mfile,'Z_scale='+strtrim(string((*pb_meta).Z_scale,format='(f10.3)'),2)
  printf,mfile,'X_offset='+strtrim(string((*pb_meta).X_offset,format='(f10.3)'),2)
  printf,mfile,'Y_offset='+strtrim(string((*pb_meta).Y_offset,format='(f10.3)'),2)
  printf,mfile,'Z_offset='+strtrim(string((*pb_meta).Z_offset,format='(f10.3)'),2)
  printf,mfile,'I_Scale='+strtrim(string((*pb_meta).I_Scale,format='(f10.3)'),2)

  flush,mfile
  
  ; Get path and file name as separate strings
  last=strpos(point_file,path_sep(),/reverse_search)
  in_path = file_dirname(point_file)
  in_base=strtrim(strmid(point_file,last+1,strlen(point_file)-last-1),2)
  
  (*pb_meta).point_file_Path=strtrim(in_path,2)
  (*pb_meta).point_file_Name=strtrim(in_base,2)
  
  flush,mfile
  
  printf,mfile,'Ptcl_File_Path='+strtrim((*pb_meta).point_file_Path,2)
  printf,mfile,'Ptcl_File_Name='+strtrim((*pb_meta).point_file_Name,2)
  
  ;===============================================================
  ;code to save solutions for d for debug
  if (ovdebug)then begin
    temp_name=strtrim(in_path,2)+path_sep()+'temp_file_'+strtrim(string(wavelength),2) $
      +'_'+strtrim(string(num_test),2)+'pts_ptcl_sol.txt'
    print,strtrim(temp_name,2)
    ;see if the TEMP file exists & remove if it does!
    if(file_test(temp_name)) then begin
      fids=envi_get_file_ids()
      if(fids[0] eq -1) then begin
        file_delete, temp_name,/quiet
        print,'old temp file deleted'
      endif else begin
        for i=0,n_elements(fids)-1 do begin
          envi_file_query,fids[i],fname=tname
          if (strtrim(strlowcase(temp_name),2) eq $
            strtrim(strlowcase(tname),2)) then begin
            envi_file_mng,id=fids[i],/remove
            print,'old temp file removed from ENVI'
          endif
        endfor
        file_delete, temp_name,/quiet
        print,'old temp file deleted'
      endelse
    endif
    text_err=0
    openw, tempfile, temp_name,/get_lun,error=text_err
    if (text_err ne 0) then begin
      print,'Error opening temp file in point cloud!!'
      print,'File Name =',strtrim(temp_name,2)
      print,'text_err=',text_err
      print,'Error Type =',strtrim(string(!ERROR_STATE.MSG),2)
      print,'Sys_Error Type =',strtrim(string(!ERROR_STATE.SYS_MSG),2)
      error=4
      goto, cleanup
    endif
    printf,tempfile,'File for selected outputs to monitor overlap correction'
    printf,tempfile,strtrim('Run made at: '+time_date,2)
    printf,tempfile,'num_test='+strtrim(string(num_test),2)
    printf,tempfile,'test_num='+strtrim(string(test_num),2)
    flush,tempfile
  endif
  ;
  ;=======================================================
  
  ;only use points with positive range greater than focal length
  pos_r=where((*pb_stats).range le 0.50,nzero)
  h=(*pb_stats).range[1]-(*pb_stats).range[0]
  if (h le 0.0) then h=1.0
  rmax=float((*pb_stats).range[n_elements((*pb_stats).range)-1])
  if (nzero gt 0) then rmax=float((*pb_stats).range[pos_r[nzero-1]]+h)
  
  print,'rmax=',rmax
  
  ; Calculate some values needed in the convolution
  p=float(p)
  psum=total(p*p)
  pnorm=total(p)
  pmax=max(p)
  spfac=psum/(pnorm*pmax)
  m = max(p,mpos)
  nlead = mpos
  ntrail = n_elements(p)-nlead-1
  
  ;NOTE: this is a change of one value (it was ntrail = n_elements(p)-nlead)
  ;print,'max p=',m
  ;print,'fwhm p=',pnorm/m
  ;print,'nlead,ntrail=',nlead,ntrail
  ;
  ;now put in code to ensure the result if there is truncation
  ; Pad pulse to make symmetrical about the peak
  ppad=float(p)
  if (ntrail gt nlead) then begin
    ppad = [replicate(0.0,ntrail-nlead),p]
  endif else if (nlead gt ntrail) then begin
    ppad = [p,replicate(0.0,nlead-ntrail)]
  endif
  
  value=max(ppad,max_point)
  ipscal=float(ppad[max_point])/mean(ppad[max_point-1:max_point+1])
  print,'ipscal=',ipscal
  print,'elements in ppad=',n_elements(ppad)
  print,'mid point=',n_elements(ppad)/2
  print,'max_point=',max_point
  
  num_pad=n_elements(ppad)
  fnum=float(num_pad)
  ;
  ppad2=[replicate(0,num_pad),ppad,replicate(0,num_pad)]
  cross2=convol(ppad2,ppad,pnorm)
  cross2=cross2[num_pad:2*num_pad-1]/pmax
  value2=max(cross2,max_point_2)
  
  print,'num_pad='+strtrim(string(num_pad),2)
  print,'size of cross2='+strtrim(string(n_elements(cross2)),2)
  print,'max in ppad=',value
  print,'max in cross2=',value2
  print,'max point in ppad=',max_point
  print,'max point in cross2=',max_point_2
  
  ;=====================================
  ;find where the left and right 0.0001*peak is
  tmp_I = where(cross2/value2 ge 0.0001)
  ; find the trough and secondary peak
  ; find the three zero-cross points after the maximum peak
  wflen=n_elements(cross2)
  zero_xloc = where(cross2[max_point_2:wflen-2]*cross2[max_point_2+1:wflen-1] le 0, tmpcount) + max_point_2
  ; find the minimum and maximum between the first zero-cross point and the third zero-cross point
  tmp = min(cross2[zero_xloc[0]:zero_xloc[2]],tmploc)
  troughloc = round(mean(tmploc)) + zero_xloc[0]
  tmp = max(cross2[zero_xloc[0]:zero_xloc[2]],tmploc)
  scdpeakloc = round(mean(tmploc)) + zero_xloc[0]
  i_val=[0,tmp_I[0],max_point_2,troughloc,scdpeakloc,tmp_I[n_elements(tmp_I)-1],n_elements(cross2)-1]
  ;=====================================
  
  print,'critical points of cross2 (i_val):'
  print,'['+strtrim(strjoin(strtrim(string(i_val),2),',',/single),2)+']'
  
  print,'critical values of cross2:'
  print,'['+strtrim(strjoin(strtrim(string(cross2[i_val]),2),',',/single),2)+']'
  
  h1=i_val[3]-i_val[2]
  h2=i_val[4]-i_val[2]
  
  print,'new h1,h2='+'['+strtrim(string(h1),2)+','+strtrim(string(h2),2)+']'
  ;save the pulse information
  ord=(findgen(num_pad)-float(max_point))*h
  min_ord=min(ord)
  max_ord=max(ord)
  value0=float(ppad)/pnorm
  value1=float(ppad)/pmax
  value2=float(cross2)
  d_fwhm=1.0/max(value0)
  p_fwhm=d_fwhm*h
  cross2=0b
  l_file=strtrim((*pb_stats).outfile,2)
  if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
    pp_file=strtrim(l_file,2)+'_pulse.txt'
  endif else begin
    pp_file=strtrim(strmid(l_file,0,n_dot),2)+'_pulse.txt'
  endelse
  
  ;see if the pp file exists & remove if it does!
  if(file_test(pp_file)) then begin
    fids=envi_get_file_ids()
    if(fids[0] eq -1) then begin
      file_delete, pp_file,/quiet
      print,'old pulse file deleted'
    endif else begin
      for i=0,n_elements(fids)-1 do begin
        envi_file_query,fids[i],fname=tname
        if (strtrim(strlowcase(pp_file),2) eq $
          strtrim(strlowcase(tname),2)) then begin
          envi_file_mng,id=fids[i],/remove
          print,'old pulse file removed from ENVI'
        endif
      endfor
      file_delete, pp_file,/quiet
      print,'old pulse file deleted'
    endelse
  endif
  text_err=0
  openw, ppfile, pp_file,/get_lun,error=text_err
  if (text_err ne 0) then begin
    print,'Error opening pulse file in point cloud!!'
    print,'File Name =',strtrim(pp_file,2)
    print,'text_err=',text_err
    print,'Error Type =',strtrim(string(!ERROR_STATE.MSG),2)
    print,'Sys_Error Type =',strtrim(string(!ERROR_STATE.SYS_MSG),2)
    error=5
    goto, cleanup
  endif
  
  ;if saving B and R images get names!
  if (save_br) then begin
    if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
      b_file=strtrim(l_file,2)+'_B_image.img'
      r_file=strtrim(l_file,2)+'_R_image.img'
    endif else begin
      b_file=strtrim(strmid(l_file,0,n_dot),2)+'_B_image.img'
      r_file=strtrim(strmid(l_file,0,n_dot),2)+'_R_image.img'
    endelse
    ;see if the b file exists & remove if it does!
    if(file_test(b_file)) then begin
      fids=envi_get_file_ids()
      if(fids[0] eq -1) then begin
        file_delete, b_file,/quiet
        print,'old B file deleted'
      endif else begin
        for i=0,n_elements(fids)-1 do begin
          envi_file_query,fids[i],fname=tname
          if (strtrim(strlowcase(b_file),2) eq $
            strtrim(strlowcase(tname),2)) then begin
            envi_file_mng,id=fids[i],/remove
            print,'old B file removed from ENVI'
          endif
        endfor
        file_delete, b_file,/quiet
        print,'old B file deleted'
      endelse
    endif
    text_err=0
    openw, bfile, b_file,/get_lun,error=text_err
    if (text_err ne 0) then begin
      print,'Error opening B file in point cloud!!'
      print,'File Name =',strtrim(b_file,2)
      print,'text_err=',text_err
      print,'Error Type =',strtrim(string(!ERROR_STATE.MSG),2)
      print,'Sys_Error Type =',strtrim(string(!ERROR_STATE.SYS_MSG),2)
      error=6
      goto, cleanup
    endif
    ;see if the R file exists & remove if it does!
    if(file_test(r_file)) then begin
      fids=envi_get_file_ids()
      if(fids[0] eq -1) then begin
        file_delete, r_file,/quiet
        print,'old R file deleted'
      endif else begin
        for i=0,n_elements(fids)-1 do begin
          envi_file_query,fids[i],fname=tname
          if (strtrim(strlowcase(r_file),2) eq $
            strtrim(strlowcase(tname),2)) then begin
            envi_file_mng,id=fids[i],/remove
            print,'old R file removed from ENVI'
          endif
        endfor
        file_delete, r_file,/quiet
        print,'old R file deleted'
      endelse
    endif
    text_err=0
    openw, rfile, r_file,/get_lun,error=text_err
    if (text_err ne 0) then begin
      print,'Error opening R file in point cloud!!'
      print,'File Name =',strtrim(r_file,2)
      print,'text_err=',text_err
      print,'Error Type =',strtrim(string(!ERROR_STATE.MSG),2)
      print,'Sys_Error Type =',strtrim(string(!ERROR_STATE.SYS_MSG),2)
      error=7
      goto, cleanup
    endif
  endif
  
  ;if saving pfilter image get name!
  if (save_pfilt) then begin
    if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
      pfilt_file=strtrim(l_file,2)+'_pfilter.img'
    endif else begin
      pfilt_file=strtrim(strmid(l_file,0,n_dot),2)+'_pfilter.img'
    endelse
    ;see if the pfilt file exists & remove if it does!
    if(file_test(pfilt_file)) then begin
      fids=envi_get_file_ids()
      if(fids[0] eq -1) then begin
        file_delete, pfilt_file,/quiet
        print,'old pfilter file deleted'
      endif else begin
        for i=0,n_elements(fids)-1 do begin
          envi_file_query,fids[i],fname=tname
          if (strtrim(strlowcase(pfilt_file),2) eq $
            strtrim(strlowcase(tname),2)) then begin
            envi_file_mng,id=fids[i],/remove
            print,'old pfilter file removed from ENVI'
          endif
        endfor
        file_delete, pfilt_file,/quiet
        print,'old pfilter file deleted'
      endelse
    endif
    text_err=0
    openw, pfile, pfilt_file,/get_lun,error=text_err
    if (text_err ne 0) then begin
      print,'Error opening pfilter file in point cloud!!'
      print,'File Name =',strtrim(pfilt_file,2)
      print,'text_err=',text_err
      print,'Error Type =',strtrim(string(!ERROR_STATE.MSG),2)
      print,'Sys_Error Type =',strtrim(string(!ERROR_STATE.SYS_MSG),2)
      error=66
      goto, cleanup
    endif
  endif
  
  ;write out information to the point cloud file - 3 headers
  time_date=strtrim(systime(),2)
  printf,ppfile,strtrim('Point Cloud Pulse Filter Information File',2)
  printf,ppfile,strtrim('Run made at: '+time_date,2)
  printf,ppfile,strtrim('num_pad=',2)+strtrim(string(num_pad),2)
  ;; fwhm from iterated pulse model in unit of samples
  printf,ppfile,strtrim('d_fwhm(samp)=',2)+strtrim(string(d_fwhm),2)
  ;; fwhm from iterated pulse model in unit of meters
  printf,ppfile,strtrim('p_fwhm(m)=',2)+strtrim(string(p_fwhm),2)
  printf,ppfile,strtrim('num,h,Pulse,I,I2',2)
  for k=0,num_pad-1 do begin
    buf=string(k+1)+','+string(ord[k])+','+string(value0[k])+','+string(value1[k])+','+string(value2[k])
    buf=strcompress(buf,/remove_all)
    printf,ppfile,strtrim(buf,2)
  endfor
  free_lun,ppfile,/force
  
  ;open the image file to save info
  ;Open output file
  text_err=0
  openw, oifile, (*pb_stats).oi_name,/get_lun,error=text_err
  if (text_err ne 0) then begin
    print,'Error opening output file '+strtrim((*pb_stats).oi_name,2)
    print,'Local error number='+strtrim(string(text_err),2)
    error=8
    goto, cleanup
  endif
  
  ;set up accumulators for the ancillary image file
  oi_accum=lonarr(nsamples,nlines)
  gaps=bytarr(nsamples,nlines)
  sum_accum=fltarr(nsamples,nlines)
  range_mean=fltarr(nsamples,nlines)
  i2_accum=fltarr(nsamples,nlines)
  d0_accum=fltarr(nsamples,nlines)
  d_accum=fltarr(nsamples,nlines)
  resid=fltarr(nsamples,nlines)
  first_hit=fltarr(nsamples,nlines)
  last_hit=fltarr(nsamples,nlines)
  mean_z=fltarr(nsamples,nlines)
  
  ;normalisation factor for correlation
  psum=sqrt(psum/fnum)
  
  ;set up and go over the image extracting the point cloud
  
  Zero_Hit_Number=0L
  Shot_Hit_Number=0L
  Total_Hit_Number=0L
  Min_Intensity=1.0e7
  Max_Intensity=-1.0e7
  Min_X=1.0e7
  Max_X=-1.0e7
  Min_Y=1.0e7
  Max_Y=-1.0e7
  Min_Z=1.0e7
  Max_Z=-1.0e7
  Shot_Num=0L
  
  num_spec=0L
  
  ;check if you are going to add the DWEL position records
  ;these are just the top and bottom of the DWEL so they are they only for geometry
  if ((*pb_meta).Add_DWEL gt 0) then begin
    buf=string(0.0,0.0,0.0,0.0,0,0,0,DWEL_num,0.0,0.0,0.0,0.0,0,0,0,format='(3f14.3,f14.4,2i10,2i14,4f14.3,3i10)')
    buf=strtrim(strcompress(buf),2)
    while (((ii = strpos(buf, ' '))) ne -1) do $
      strput, buf, ',', ii
    printf,tfile,buf
    buf=string(0.0,0.0,(*pb_meta).DWEL_Height,0.0,0,0,0,DWEL_num,0.0,0.0,0.0,0.0,0,0,0,format='(3f14.3,f14.4,2i10,2i14,4f14.3,3i10)')
    buf=strtrim(strcompress(buf),2)
    while (((ii = strpos(buf, ' '))) ne -1) do $
      strput, buf, ',', ii
    printf,tfile,buf
    flush,tfile
    Total_Hit_Number=2L
  endif
  
  accumn=0.0d0
  accumd=0.0d0
  iprint=0L
  
  ;open the input file
  err=0
  ;open input file
  openr,inlun,infile,/get_lun,error=err
  if (err gt 0) then begin
    print,'Error opening input file'
    print,'File Name='+strtrim(infile,2)
    print,'Local error='+strtrim(string(err),2)
    error=9
    goto,cleanup
  endif

  pulse = (*pb_stats).pulse
  p_range = (*pb_stats).p_range
  pulse_len = n_elements(pulse)
  ;set up the pointer for read_binary
  bufrs=long64(nbytes)*long64(nbands)*long64(nsamples)
  pointsz=long64(0)
  ;
  bad_peaks=0L
  bad_interpol=0L
  nump=0
  nump_new=0
  ; Loop through each line, calculate B and apply filter.
  ; Write to output file.
  for i=0,nlines-1 do begin 
    ;    inline = envi_get_slice(fid=fid, /bil, pos=pos, line=i, xs=0, xe=nsamples-1)
    pointsz=long64(i)*long64(bufrs)
    inline=read_binary(inlun,data_start=pointsz,data_dims=[nsamples,nbands],data_type=dt)
    if (save_br) then begin
      B_Mat=fltarr(nsamples,nbands)
      R_Mat=fltarr(nsamples,nbands)
    endif
    if (save_pfilt) then begin
      P_Mat=make_array(nsamples,nbands,type=dt)
    endif
    ;now go over the BIL slice for each sample and find points
   for j=0, nsamples-1 do begin 
      b=fltarr(nbands)
      db=fltarr(nbands)
      d2b=fltarr(nbands)
      r=fltarr(nbands)
      dr=fltarr(nbands)
      d2r=fltarr(nbands)
      if ((*pb_stats).mask[j,i] le 0b) then begin
        oi_accum[j,i]=0l
        sum_accum[j,i]=0.0
        range_mean[j,i]=0.0
        i2_accum[j,i]=0.0
        d0_accum[j,i]=0.0
        d_accum[j,i]=0.0
        resid[j,i]=0.0
        first_hit[j,i]=0.0
        last_hit[j,i]=0.0
        gaps[j,i]=0b
        mean_z[j,i]=0.0
        goto,skip
      endif
      shot_num=shot_num+1L
      ;Pad data to allow convolution to the ends of the original array
      ;note using mean in padding areas
      temp=float(reform(inline[j,*]))
      t=[replicate(mean(temp[0:19]),num_pad),temp,replicate(mean(temp[nbands-20:nbands-1]),num_pad)]
      temp=0b
      c = convol(t,ppad,pnorm)
      dc=deriv(c)
      d2c=deriv(dc)
      b = c[num_pad:num_pad+nbands-1]
      db=dc[num_pad:num_pad+nbands-1]/h
      d2b = d2c[num_pad:num_pad+nbands-1]/(h^2)
      ;
      ;get moving average value variance to compute correlation
      temp = sqrt(float(smooth(t^2,num_pad)))
      temp = temp[num_pad:num_pad+nbands-1]
      t = t[num_pad:num_pad+nbands-1]
      ;
      bs1 = shift(temp,1)
      bsm1 = shift(temp,-1)
      test=abs(temp) le b_thresh
      test1=abs(bs1) le b_thresh
      test2=abs(bsm1) le b_thresh
      w = where((test and test1 and test2), nw, compl=fok, ncompl=nfok)
      if (nw gt 0) then begin
        r[w] = 0.0
        dr[w]=0.0
        d2r[w]=0.0
        if (nfok gt 0) then begin
          r[fok] = b[fok]*psum/temp[fok]
          dr[fok] = db[fok]*psum/temp[fok]
          d2r[fok] = d2b[fok]*psum/temp[fok]
        endif
      endif else begin
        r = b*psum/temp
        dr=db*psum/temp
        d2r=d2b*psum/temp
      endelse
      ;temp=0b
      ; Check neighbourhoods of derivative and second derivative of correlation for peaks
      bs1 = shift(dr,1)
      bs2=shift(dr,2)
      bsm1 = shift(dr,-1)
      bsm2=shift(dr,-2)
      test1 = (bs1 gt 0.0001 and bs2 gt 0.0001)
      testm1 = (bsm1 lt -0.0001 and bsm2 lt -0.0001)
      test = (d2r lt -0.01) and (r gt 0.1)
      peaks = where(((test) and (test1 and testm1)),nump)
      nump_new=nump
      ph=(*pb_stats).azimuth[j,i]
      th=(*pb_stats).zenith[j,i]
      offset=0s
      ;check for runs of peaks and simplify (can happen)
      ;then score the peaks and eliminate likely non-peaks
      if (nump gt 0) then begin
        if (nump gt 1) then begin
          for k=0,nump-2 do begin
            nn=nump-k-1
            if (peaks[nn-1] eq peaks[nn]-1) then begin
              if(r[peaks[nn-1]] ge r[peaks[nn]]) then begin
                peaks[nn]=peaks[nn-1]
              endif
              peaks[nn-1:nump-2]=peaks[nn:nump-1]
              nump_new=nump_new-1
            endif
          endfor
          peaks=peaks[0:nump_new-1]
        endif
        if (nump_new gt 0) then begin
          kk=0
          peaks_out=peaks
          if (nump_new eq 1) then small_thresh=cvthresh*b[peaks[0]] else small_thresh=cvthresh*max(b[peaks])
          small_thresh=max([small_thresh,sieve_thresh])
          for k=0,nump_new-1 do begin
            score=0
            rkv=(*pb_stats).range[peaks[k]]
            zpeak=rkv*cos(th*!dtor)+(*pb_meta).DWEL_Height
            xpeak=rkv*sin(th*!dtor)*sin(ph*!dtor)
            ypeak=rkv*sin(th*!dtor)*cos(ph*!dtor)
            ;
            if ((b[peaks[k]] gt small_thresh) and (r[peaks[k]] gt r_thresh) $
              and (zpeak lt zhigh) and (zpeak gt zlow) $
              and (xpeak lt xmax) and (xpeak gt xmin) $
              and (ypeak lt ymax) and (ypeak gt ymin) $
              and (rkv gt rmax)) then begin
              posp1=min([peaks[k]+h1,nbands-1])
              posp2=min([peaks[k]+h2,nbands-1])
              posm1=max([peaks[k]-h1,0])
              posm2=max([peaks[k]-h2,0])
              if (b[posp1] lt 0.0) then score=score+1
              if (b[posm1] lt 0.0) then score=score+1
              if (b[posp2] gt 0.0) then begin
                score=score+1
                if (b[posp2] lt b[peaks[k]]) then score=score+2
              endif
              if (b[posm2] gt 0.0) then begin
                score=score+1
                if (b[posm2] lt b[peaks[k]]) then score=score+2
              endif
            endif
            if (score ge 6) then begin
              peaks_out[kk]=peaks[k]
              kk=kk+1
            endif
          endfor
          if (kk gt 0) then peaks=reform(peaks_out[0:kk-1]) else peaks=0b
          nump_new=kk
        endif
        if (nump_new le 0) then goto,nohits
        ;
        ;Now look over the selected peaks and get (x,y,z,I) values
        ;
        redone=0b
        solveit:
        esum=0.0
        bsum=0.0
        rg=fltarr(nump_new)
        intensity=fltarr(nump_new)
        I2=fltarr(nump_new)
        speaks=fltarr(nump_new)
        for k=0,nump_new-1 do begin
          val=r[peaks[k]-1:peaks[k]+1]
          if (val[1] lt (val[0]+val[2])/2.0) then begin
            bad_peaks=long(bad_peaks)+1L
          endif
          rg_loc=(*pb_stats).range[peaks[k]-1:peaks[k]+1]
          istat=peak_int(rg_loc,val,rg_peak,intensity_peak,offset)
          if (istat gt 0) then bad_interpol=long(bad_interpol)+1L
          rg[k]=rg_peak
          speaks[k]=peaks[k]+offset
          intensity[k]=ipscal*mean(t[peaks[k]-1:peaks[k]+1])
          I2[k]=b[speaks[k]]
          esum=esum+float(intensity[k])
          bsum=bsum+float(I2[k])
          if ((nump_new eq 1) and (abs(I2[k]) gt 0.0001)) then begin
            accumn=accumn+double(I2[k])^2
            accumd=accumd+double(intensity[k])^2
          endif
        endfor
        ;
        ;Next correct for overlap and use I and I2 to estimate stable values
        ls_res=0.0
        if (nump_new gt 0) then begin
          difmat=cmreplicate(rg,nump_new)
          difmat=difmat-transpose(difmat)
          pos_out=where((difmat lt min_ord) or (difmat gt max_ord),num_out)
          ;here interpolate I and I2 with dismat, set zeros and then you are close.
          ;two matrices are done and vectors are I and I2 then multiply
          I_mat=interpol(value1,ord,difmat)
          I2_mat=interpol(value2,ord,difmat)
          if (num_out gt 0) then begin
            I_mat[pos_out]=0.0
            I2_mat[pos_out]=0.0
          endif
          I_mat=reform(I_mat,nump_new,nump_new)
          I2_mat=reform(I2_mat,nump_new,nump_new)
          ls_mat=fltarr(nump_new,2*nump_new)
          ls_rhs=fltarr(2*nump_new)
          ls_mat=[[I_mat],[I2_mat]]
          ls_rhs=[intensity,I2]
          d_out=la_least_squares(ls_mat,ls_rhs,residual=ls_res)
          ls_res=sqrt(abs(float(ls_res))/float(2*nump_new))
          
          if (~redone) then begin
            redone=1b
            posr=where(d_out gt small_thresh,nposr)
            if (nposr le 0) then begin
              nump_new=0
              peaks=0b
              goto,nohits
            endif
            if (nposr lt nump_new) then begin
              peaks=reform(peaks[posr])
              nump_new=nposr
              goto,solveit
            endif
            posr=0b
          endif
          ;
          ;Code to save information for debugging the linear solution for overlap
          ;======================================================================
          saving=0b
          if (ovdebug) then begin
            if ((nump_new ge num_test) and (iprint le test_num) and (ls_res gt $
              3.5)) then begin
;;              and (j gt 100) and (j lt 200) and (i gt 50) $
;;              ) then begin
              printf,tempfile,''
              printf,tempfile,'Test Case Number='+strtrim(string(iprint+1),2)
              printf,tempfile,'Sample='+strtrim(string(j+1),2)
              printf,tempfile,'Line='+strtrim(string(i+1),2)
              flush,tempfile
              printf,tempfile,'rg range vector'
              buf=strtrim(string(rg),2)
              buf=strcompress(buf,/remove_all)
              buf=strtrim(strjoin(buf,' '),2)
              printf,tempfile,buf
              buf=''
              printf,tempfile,'difmat'
              for jj=0,nump_new-1 do begin
                buf=strtrim(string(reform(difmat[*,jj])),2)
                buf=strcompress(buf,/remove_all)
                buf=strtrim(strjoin(buf,' '),2)
                printf,tempfile,buf
                buf=''
              endfor
              flush,tempfile
              printf,tempfile,'I_mat'
              for jj=0,nump_new-1 do begin
                buf=strtrim(string(reform(I_mat[*,jj])),2)
                buf=strcompress(buf,/remove_all)
                buf=strtrim(strjoin(buf,' '),2)
                printf,tempfile,buf
                buf=''
              endfor
              flush,tempfile
              printf,tempfile,'I2_mat'
              for jj=0,nump_new-1 do begin
                buf=strtrim(string(reform(I2_mat[*,jj])),2)
                buf=strcompress(buf,/remove_all)
                buf=strtrim(strjoin(buf,' '),2)
                printf,tempfile,buf
                buf=''
              endfor
              flush,tempfile
              printf,tempfile,'ls_mat'
              for jj=0,2*nump_new-1 do begin
                buf=strtrim(string(reform(ls_mat[*,jj])),2)
                buf=strcompress(buf,/remove_all)
                buf=strtrim(strjoin(buf,' '),2)
                printf,tempfile,buf
                buf=''
              endfor
              flush,tempfile
              printf,tempfile,'ls_rhs'
              buf=strtrim(string(ls_rhs),2)
              buf=strcompress(buf,/remove_all)
              buf=strtrim(strjoin(buf,' '),2)
              printf,tempfile,buf
              buf=''
              flush,tempfile
              printf,tempfile,'d_out'
              buf=strtrim(string(d_out),2)
              buf=strcompress(buf,/remove_all)
              buf=strtrim(strjoin(buf,' '),2)
              printf,tempfile,buf
              buf=''
              flush,tempfile
              printf,tempfile,'ls_res'
              buf=strtrim(string(ls_res),2)
              buf=strtrim(strjoin(buf,' '),2)
              printf,tempfile,buf
              buf=''
              flush,tempfile
              if (~(*pb_stats).cal_dat) then printf,tempfile,''
              iprint=iprint+1L
              saving=1b
            endif
          endif
          ;
          ;======================================================================
          
          difmat=0b
          I_mat=0b
          I2_mat=0b
          ls_mat=0b
          ls_rhs=0b
        endif
        if (nump_new gt noise_hits) then begin
          nump_new=0
          goto,nohits
        endif
        d0_out=reform(d_out[0:nump_new-1])
        d0sum=0.0
        for k=0,nump_new-1 do begin
          d0sum=d0sum+float(d0_out[k])
        endfor
        rg=reform(rg[0:nump_new-1])
        ;
        ;Now calibrate the d values and get the ls_res (residual of LS) if required
        if ((*pb_stats).cal_dat and ~DWEL_AppRefl) then begin
          if strcmp(laser_man, 'manlight') then begin
            eff=DWEL_eff_nsf(wavelength,rg, par=(*pb_stats).eff_par)
          endif
          if strcmp(laser_man, 'keopsys') then begin
            eff=DWEL_eff_oz(wavelength,rg, par=(*pb_stats).eff_par)
          endif

          temp=(*pb_stats).s_Factor*(rg^(*pb_stats).rpow)*d_out/(eff*(*pb_stats).DWEL_cal)
          if ((n_elements(temp) ne n_elements(d_out)) or $
              (n_elements(temp) ne n_elements(rg))) then begin
            print,'Bad error! temp and d_out do NOT conform!'
            print,'number of elements in temp=',n_elements(temp)
            print,'number of elements in d_out=',n_elements(d_out)
            print,'number of elements in rg=',n_elements(rg)
          endif
          d_out=temp
          temp=0b
          temp=(*pb_stats).s_Factor*(rg^(*pb_stats).rpow)*intensity/(eff*(*pb_stats).DWEL_cal)
          intensity = temp
          temp = 0b

          if strcmp(laser_man, 'manlight') then begin
            eff=DWEL_eff_nsf(wavelength,(*pb_stats).range, par=(*pb_stats).eff_par)
          endif
          if strcmp(laser_man, 'keopsys') then begin
            eff=DWEL_eff_oz(wavelength,(*pb_stats).range, par=(*pb_stats).eff_par)
          endif

          rvalid = where((*pb_stats).range gt rmax, nvalid)
          temp = fltarr(size(inline[j,*], /n_elements))
          if nvalid gt 0 then begin
            tmprange = (*pb_stats).range
            temp[rvalid] = $
            (*pb_stats).s_Factor*tmprange[rvalid]^(*pb_stats).rpow*reform(inline[j, $
              rvalid])/(eff[rvalid]*(*pb_stats).DWEL_cal)
            inline[j, *] = fix(round(temp*i_scale), type=2)
          endif 
          temp = 0b
          if (ovdebug and saving) then begin
            printf,tempfile,''
            printf,tempfile,'eff'
            buf=strtrim(string(eff),2)
            buf=strtrim(strjoin(buf,' '),2)
            printf,tempfile,buf
            buf=''
            flush,tempfile 
;            printf,tempfile,''
            printf,tempfile,'d_out'
            buf=strtrim(string(d_out),2)
            buf=strtrim(strjoin(buf,' '),2)
            printf,tempfile,buf
            buf=''
            flush,tempfile
            printf,tempfile,''
          endif
        endif
        ;
        rsum=0.0
        vsum=0.0
        r0sum=0.0
        v0sum=0.0
        for k=0,nump_new-1 do begin
          rsum=rsum+rg[k]*float(d_out[k])
          vsum=vsum+float(d_out[k])
          r0sum=r0sum+rg[k]*float(d0_out[k])
          v0sum=v0sum+float(d0_out[k])
        endfor
        ;
        if ((*pb_stats).cal_dat and ((vsum le -0.01) or (vsum ge 2.5))) then begin
          oi_accum[j,i]=0l
          sum_accum[j,i]=0.0
          range_mean[j,i]=0.0
          i2_accum[j,i]=0.0
          d0_accum[j,i]=0.0
          d_accum[j,i]=0.0
          resid[j,i]=0.0
          first_hit[j,i]=0.0
          last_hit[j,i]=0.0
          gaps[j,i]=0b
          mean_z[j,i]=0.0
          (*pb_stats).mask[j,i]=0b
          num_spec=num_spec+1L
          goto,skip
        endif else begin
          Shot_Hit_Number=Shot_Hit_Number+1L
          Total_Hit_Number=Total_Hit_Number+nump_new
          oi_accum[j,i]=long(nump_new)
          sum_accum[j,i]=float(esum)
          i2_accum[j,i]=float(bsum)
          d0_accum[j,i]=float(d0sum)
          d_accum[j,i]=float(vsum)
          range_mean[j,i]=r0sum/v0sum
          first_hit[j,i]=float(rg[0])
          last_hit[j,i]=float(rg[nump_new-1])
          gaps[j,i]=0b
          resid[j,i]=ls_res
          mean_z[j,i]=range_mean[j,i]*cos(th*!dtor)+(*pb_meta).DWEL_Height
        endelse
        
        return_fwhm = 0b
        return_fwhm = fltarr(nump_new)
        ; now calculate integral of return pulse of each point and the FWHM
        for k=0,nump_new-1 do begin
          ;; to see if there is another peak in close range within the range
          ;; extent of integral          
          range_left=rg[k]-range_extent*2
          range_right=rg[k]+range_extent*2
          cpind = where(rg gt range_left and rg lt range_right, ncp)
          if ncp gt 0 then tmpind = where(cpind ne k, ncp)
          ;;inline is scaled by i_scale, need to de-scale here.
          tmpwf = inline[j, *]/float(i_scale)
          if ncp gt 0 then begin
            cpind = cpind[tmpind]
            for icp = 0, ncp-1 do begin
              ;; subtract each close return pulse from the waveform before
              ;; calculating integral of this peak
              in_range = p_range + rg[cpind[icp]]
              out_range_pos = indgen(pulse_len) + (peaks[cpind[icp]]-fix(pulse_len/2))
              out_range = (*pb_stats).range[out_range_pos]
              cp_out = interpol(d_out[cpind[icp]]*pulse, in_range, out_range)
              tmpwf[out_range_pos] = tmpwf[out_range_pos] - cp_out
            endfor
          endif 
          range_left=rg[k]-range_extent
          range_right=rg[k]+range_extent
          posext=where(((*pb_stats).range ge range_left) and ((*pb_stats).range le range_right),nposext)
          if nposext gt 0 then begin
            return_fwhm[k] = total(tmpwf[posext])/d_out[k]
          endif 
        endfor 

        ;Now go over the final point cloud model and write out the records to the point cloud file
        for k=0,nump_new-1 do begin
          x=rg[k]*sin(th*!dtor)*sin(ph*!dtor)
          y=rg[k]*sin(th*!dtor)*cos(ph*!dtor)
          z=rg[k]*cos(th*!dtor)+(*pb_meta).DWEL_Height
          buf=string(x,y,z,i_scale*d_out[k],k+1,nump_new,shot_num,DWEL_num,rg[k],th,ph,(*pb_stats).range[peaks[k]],j+1,i+1,peaks[k]+1,return_fwhm[k],format='(3f14.3,f14.4,2i10,2i14,4f14.3,3i10,f14.3)')
          buf=strtrim(strcompress(buf),2)
          while (((ii = strpos(buf, ' '))) ne -1) do $
            strput, buf, ',', ii
          printf,tfile,buf
          if (i_scale*d_out[k] gt max_intensity) then max_intensity=i_scale*d_out[k]
          if (i_scale*d_out[k] lt min_intensity) then min_intensity=i_scale*d_out[k]
          if (x gt max_x) then max_x=x
          if (x lt min_x) then min_x=x
          if (y gt max_y) then max_y=y
          if (y lt min_y) then min_y=y
          if (z gt max_z) then max_z=z
          if (z lt min_z) then min_z=z
        endfor
      endif else begin
        ; Otherwise it is a gap
      
        nohits:
        Zero_Hit_Number=Zero_Hit_Number+1L
        oi_accum[j,i]=0l
        sum_accum[j,i]=0.0
        i2_accum[j,i]=0.0
        d0_accum[j,i]=0.0
        d_accum[j,i]=0.0
        range_mean[j,i]=0.0
        first_hit[j,i]=0.0
        last_hit[j,i]=0.0
        gaps[j,i]=1b
        resid[j,i]=0.0
        mean_z[j,i]=0.0
        if ((*pb_meta).zero_hit_option gt 0) then begin
          buf=string(0.0,0.0,0.0,0.0,0,0,shot_num,DWEL_num,0.0,th,ph,0.0,j+1,i+1,0,0,format='(3f14.3,f14.4,2i10,2i14,4f14.3,3i10,f14.3)')
          buf=strtrim(strcompress(buf),2)
          while (((ii = strpos(buf, ' '))) ne -1) do $
            strput, buf, ',', ii
          printf,tfile,buf
        endif
      endelse
      skip:
      if (save_br) then begin
        B_Mat[j,*]=b
        R_Mat[j,*]=r
      endif
      if (save_pfilt) then begin
        if (nump_new gt 0) then begin
          filtloc=make_array(nbands,type=dt)
          for k=0,nump_new-1 do begin
            range_left=(*pb_stats).range[peaks[k]]-range_extent
            range_right=(*pb_stats).range[peaks[k]]+range_extent
            posext=where(((*pb_stats).range ge range_left) and ((*pb_stats).range le range_right),nposext)
            if (nposext gt 0) then filtloc[posext]=1
          endfor
          P_Mat[j,*]=reform(inline[j,*])*filtloc
        endif
      endif
      b=0b
      db=0b
      d2b=0b
      r=0b
      dr=0b
      d2r=0b
      c=0b
      dc=0b
      d2c=0b
      bs1=0b
      bsm1=0b
      bs2=0b
      bsm2=0b
    endfor
    inline=0b
    if (save_br) then begin
      writeu,bfile,B_Mat
      writeu,rfile,R_Mat
      B_Mat=0b
      R_Mat=0b
    endif
    if (save_pfilt) then begin
      writeu,pfile,P_Mat
      P_Mat=0b
    endif
  endfor
  flush,tfile
  if (bad_peaks gt 0L) then print,'WARNING - There were '+strtrim(string(bad_peaks),2)+' bad peaks'
  if (bad_interpol gt 0L) then print,'WARNING - There were '+strtrim(string(bad_interpol),2)+' bad interpolations'
  free_lun,inlun,/force
  ;
  if (save_br) then begin
    free_lun,bfile,/force
    free_lun,rfile,/force
    b_descrip='B_image for '+strtrim(descrip,2)
    envi_setup_head,fname=b_file,ns=nsamples,nl=nlines,nb=nbands,$
      xstart=0,ystart=0,$
      data_type=4, interleave=1, $
      descrip=b_descrip, wl=wl, bnames=bnames,/write
      envi_open_file,b_file,r_fid=fid,/no_interactive_query,/no_realize
      ;write out the previous header records
      status=dwel_put_headers(fid,DWEL_headers)
      envi_file_mng,id=fid,/remove
;
    r_descrip='R_image for '+strtrim(descrip,2)
    envi_setup_head,fname=r_file,ns=nsamples,nl=nlines,nb=nbands,$
      xstart=0,ystart=0,$
      data_type=4, interleave=1, $
      descrip=r_descrip, wl=wl, bnames=bnames,/write
    envi_open_file,r_file,r_fid=fid,/no_interactive_query,/no_realize
    ;write out the previous header records
    status=dwel_put_headers(fid,DWEL_headers)
    envi_file_mng,id=fid,/remove
;
  endif
  if (save_pfilt) then begin
    free_lun,pfile,/force
    p_descrip='Pfilter_image for '+strtrim(descrip,2)
    envi_setup_head,fname=pfilt_file,ns=nsamples,nl=nlines,nb=nbands,$
      xstart=0,ystart=0,$
      data_type=dt, interleave=1, $
      descrip=p_descrip, wl=wl, bnames=bnames,/write
      envi_open_file,pfilt_file,r_fid=fid,/no_interactive_query,/no_realize
      ;write out the previous header records
      status=dwel_put_headers(fid,DWEL_headers)
      envi_file_mng,id=fid,/remove
;

  endif
  if (ovdebug) then begin
    if (iprint le 0) then begin
      printf,tempfile,'There were NO such cases detected!'
      flush,tempfile
    endif
    flush,tempfile
    free_lun,tempfile,/force
  endif
  
  print,'num_spec=',num_spec
  
  ;compute the actual ratio to compare with theory
  ratio_th=spfac
  ratio_est=sqrt(double(accumn)/double(accumd))
  pos_mask=where((*pb_stats).mask ne 0b,npos)
  tmask=(*pb_stats).mask
  
  ;get d-stats
  d_accum = d_accum * i_scale
  image_statistics, d_accum,mask=tmask,minimum=emin,maximum=emax,mean=emean,stddev=esdev
  emin = emin / float(i_scale)
  emax = emax / float(i_scale)
  emean = emean / float(i_scale)
  esdev = esdev / float(i_scale)
;;  d_accum[pos_mask]=4095.0*(d_accum[pos_mask]-emin)/(emax-emin)
  print,'mean corrected intensity=',emean
  
  pb_info=[$
    'Stats_Format=(Min,Mean,Max,Stddev)',$
    'd_Stats=('+strtrim(string(emin),2)+',' $
    +strtrim(string(emean),2)+',' $
    +strtrim(string(emax),2)+',' $
    +strtrim(string(esdev),2)+')' $
    ]
    
  ;get d-stats
  image_statistics,d0_accum,mask=tmask,minimum=emin,maximum=emax,mean=emean,stddev=esdev
;;  d0_accum[pos_mask]=4095.0*(d0_accum[pos_mask]-emin)/(emax-emin)
  print,'mean corrected uncalibrated intensity=',emean
  
  pb_info=[pb_info,$
    'd0_Stats=('+strtrim(string(emin),2)+',' $
    +strtrim(string(emean),2)+',' $
    +strtrim(string(emax),2)+',' $
    +strtrim(string(esdev),2)+')' $
    ]
    
  ;get I stats
  image_statistics,sum_accum,mask=tmask,minimum=emin,maximum=emax,mean=emean,stddev=esdev
;;  sum_accum[pos_mask]=4095.0*(sum_accum[pos_mask]-emin)/(emax-emin)
  print,'mean I sum=',emean
  
  pb_info=[pb_info,$
    'I_Stats=('+strtrim(string(emin),2)+',' $
    +strtrim(string(emean),2)+',' $
    +strtrim(string(emax),2)+',' $
    +strtrim(string(esdev),2)+')' $
    ]
    
  ;get I2 stats
  image_statistics,i2_accum,mask=tmask,minimum=emin,maximum=emax,mean=emean,stddev=esdev
;;  i2_accum[pos_mask]=4095.0*(i2_accum[pos_mask]-emin)/(emax-emin)
  print,'mean I2 sum=',emean
  
  pb_info=[pb_info,$
    'I2_Stats=('+strtrim(string(emin),2)+',' $
    +strtrim(string(emean),2)+',' $
    +strtrim(string(emax),2)+',' $
    +strtrim(string(esdev),2)+')' $
    ]
    
  ;get Range stats
  image_statistics,range_mean,mask=tmask,minimum=emin,maximum=emax,mean=emean,stddev=esdev
  print,'mean Range=',emean
  
  pb_info=[pb_info,$
    'Range_Stats=('+strtrim(string(emin),2)+',' $
    +strtrim(string(emean),2)+',' $
    +strtrim(string(emax),2)+',' $
    +strtrim(string(esdev),2)+')' $
    ]
    
  ;get residual stats
  image_statistics,resid,mask=tmask,minimum=emin,maximum=emax,mean=emean,stddev=esdev
;;  resid[pos_mask]=4095.0*(resid[pos_mask]-emin)/(emax-emin)
  print,'mean Residual=',emean
  
  pb_info=[pb_info,$
    'Residual_Stats=('+strtrim(string(emin),2)+',' $
    +strtrim(string(emean),2)+',' $
    +strtrim(string(emax),2)+',' $
    +strtrim(string(esdev),2)+')' $
    ]
    
  ;write out the ancillary file in the at-project geometry
  writeu,oifile,long(oi_accum)
  writeu,oifile,round(d_accum)
  writeu,oifile,round(d0_accum)
  writeu,oifile,round(sum_accum)
  writeu,oifile,round(i2_accum)
  writeu,oifile,round(100.0*first_hit)
  writeu,oifile,round(100.0*range_mean)
  writeu,oifile,round(100.0*last_hit)
  writeu,oifile,round(100.0*(*pb_stats).zenith)
  writeu,oifile,round(100.0*(*pb_stats).azimuth)
  writeu,oifile,round(resid)
  writeu,oifile,long(gaps)
  writeu,oifile,long((*pb_stats).mask)
  writeu,oifile,round(100.0*mean_z)
  free_lun, oifile,/force
  pos_mask=0b
  tmask=0b
  
  oi_accum=0b
  d_accum=0b
  d0_accum=0b
  sum_accum=0b
  i2_accum=0b
  first_hit=0b
  last_hit=0b
  range_mean=0b
  
  if ((*pb_meta).zero_hit_option gt 0) then begin
    Nrecs=Total_Hit_Number+Zero_Hit_Number
  endif else begin
    Nrecs=Total_Hit_Number
  endelse
  
  (*pb_meta).Zero_Hit_Number=Zero_Hit_Number
  (*pb_meta).Shot_Hit_Number=Shot_Hit_Number
  (*pb_meta).Total_Hit_Number=Total_Hit_Number
  (*pb_meta).Total_Hit_Number=Total_Hit_Number
  (*pb_meta).Nrecs=Nrecs
  (*pb_meta).Min_X=Min_X
  (*pb_meta).Max_X=Max_X
  (*pb_meta).Min_Y=Min_Y
  (*pb_meta).Max_Y=Max_Y
  (*pb_meta).Min_Z=Min_Z
  (*pb_meta).Max_Z=Max_Z
  (*pb_meta).Min_Intensity=Min_Intensity
  (*pb_meta).Max_Intensity=Max_Intensity
  
  printf,mfile,'Zero_Hit_Number='+strtrim(string((*pb_meta).Zero_Hit_Number,format='(i10)'),2)
  printf,mfile,'Shot_Hit_Number='+strtrim(string((*pb_meta).Shot_Hit_Number,format='(i10)'),2)
  printf,mfile,'Total_Hit_Number='+strtrim(string((*pb_meta).Total_Hit_Number,format='(i10)'),2)
  printf,mfile,'Nrecs='+strtrim(string((*pb_meta).Nrecs,format='(i10)'),2)
  printf,mfile,'Max_X='+strtrim(string((*pb_meta).Max_X,format='(f10.3)'),2)
  printf,mfile,'Min_X='+strtrim(string((*pb_meta).Min_X,format='(f10.3)'),2)
  printf,mfile,'Max_Y='+strtrim(string((*pb_meta).Max_Y,format='(f10.3)'),2)
  printf,mfile,'Min_Y='+strtrim(string((*pb_meta).Min_Y,format='(f10.3)'),2)
  printf,mfile,'Max_Z='+strtrim(string((*pb_meta).Max_Z,format='(f10.3)'),2)
  printf,mfile,'Min_Z='+strtrim(string((*pb_meta).Min_Z,format='(f10.3)'),2)
  printf,mfile,'Max_Intensity='+strtrim(string((*pb_meta).Max_Intensity,format='(f10.3)'),2)
  printf,mfile,'Min_Intensity='+strtrim(string((*pb_meta).Min_Intensity,format='(f10.3)'),2)
  
  flush,mfile
  
  perc_hit=(100.0*float((*pb_meta).Shot_Hit_Number)/float((*pb_meta).Total_Hit_Number))

  pb_info=[pb_info,$
    'Rmax='+strtrim(string(rmax,format='(f10.3)'),2),$
    'Percent_hits='+strtrim(string(perc_hit,format='(f10.3)'),2),$
    'Theory_Ratio='+strtrim(string(ratio_th,format='(f10.4)'),2),$
    'Estimated_Ratio='+strtrim(string(ratio_est,format='(f10.4)'),2) $
    ]
    
  print,'Point Cloud Case Finished - clean up'
  
  cleanup:
  
  ; Close output files
  free_lun, tfile,/force
  free_lun, mfile,/force
  free_lun,tempfile,/force
  free_lun,bfile,/force
  free_lun,rfile,/force
  free_lun,inlun,/force
  free_lun,pfile,/force
  free_lun,ppfile,/force
  
  envi_file_mng,id=fid,/remove
  
  return
end

;======================================================================
;; settings is a structure variable containing all setting for point cloud
;; generation. 
pro dwel_get_point_cloud, infile, ancfile, outfile, err, Settings=settings

  compile_opt idl2
;  envi, /restore_base_save_files
;  envi_batch_init, /no_status_window

  resolve_routine, 'DWEL_GET_HEADERS', /compile_full_file, /either
  resolve_routine, 'DWEL_ITPULSE_MODEL_DUAL_NSF', /compile_full_file, /either
  resolve_routine, 'DT2NB', /compile_full_file, /either
  resolve_routine, 'DWEL_PUT_HEADERS', /compile_full_file, /either
  resolve_routine, 'CMREPLICATE', /compile_full_file, /either

  ;; get the size of input file to be processed. It will be used in later
  ;; summary of processing time. 
  procfilesize = file_info(infile)
  procfilesize = procfilesize.size
  ;; get the time now as the start of processing
  starttime = systime(1)

  ; infile is the dwel file of data
  ; ancfile is its ancillary
  ; cal_dat=0 no calibration =1 carry out calibration
  ; dwel_az_n is the dwell azimuth to North
  ; runcode is a number (integer or long integer)unique to the run
  ; err is the error flag =0 OK >0 error number indicates location or type
  
  err=0
  tfile=30
  dwel_pointcloud_info=['']
  
  ;; set up a structure to store default settings and will be updated
  ;; accordingly later. 
  ;; runcode: set to a value to distinguish runs, be default it is set to the
  ;; current Julian day*10.
  ;; add_dwel: if 1, two points (0, 0, 0) and (0, 0, dwel_height) are recorded
  ;; in generated point cloud for reference.
  ;; save_br: if 1, save images of b and r. they are really really large bc they
  ;; are image of doulbe floating values. Only save them when debugging.
  ;; save_pfilt: if 1, save pfilter image
  ;; zlow, zhigh, xmin, xmax, ymin, ymax: set a bounding box of limits for
  ;; impossible or unnecessary points useful to remove impossible points
  ;; sdevfac: how many times of standard deviation of noise to determine a
  ;; threshold to filter out noise points.
  finalsettings = { ptclsettings, $
    cal_dat:0, $
    DWEL_az_n:0, $ ; unit, deg
    runcode:round(systime(/julian)*10), $
    save_zero_hits:1, $
    add_dwel:0, $
    save_br:0, $
    save_pfilt:1, $
    zlow:-5.0, $
    zhigh:50.0, $
    xmin:-50.0, $
    xmax:50.0, $
    ymin:-50.0, $
    ymax:50.0, $
    sdevfac:2.0, $
    r_thresh:0.175, $
    sievefac:10.0, $
    cal_par:[-1,-1,-1,-1,-1,-1] $ ;[c0, c1, c2, c3, c4, b]
    }
  ;; tag names we need in settings
  setting_tag_names = tag_names(finalsettings)
  if n_elements(settings) ne 0 or arg_present(settings) then begin
    ;; user supplied settings
    ;; go through user supplied settings and update the final settings. 
    numtags = n_tags(settings)
    tags = tag_names(settings)
    for n=0, numtags-1 do begin
      tmpind = where(strmatch(setting_tag_names, tags[n], /fold_case) eq 1, $
        tmpnum) 
      if tmpnum eq 1 then begin
        finalsettings.(tmpind) = settings.(n)
      endif else begin
        print, 'Tag name is invalid or ambiguous in given ' + $
          'settings. Default value will be used instead. '
        print, 'Given tag name = ' + strtrim(tags[n], 2)
      endelse 
    endfor 
  endif 

  cal_dat=finalsettings.cal_dat
  DWEL_az_n=finalsettings.DWEL_az_n
  runcode=finalsettings.runcode
  save_zero_hits=finalsettings.save_zero_hits
  add_dwel=finalsettings.add_dwel
  save_br=finalsettings.save_br
  save_pfilt=finalsettings.save_pfilt
  zlow=finalsettings.zlow
  zhigh=finalsettings.zhigh
  xmin=finalsettings.xmin
  xmax=finalsettings.xmax
  ymin=finalsettings.ymin
  ymax=finalsettings.ymax
  sdevfac=finalsettings.sdevfac
  r_thresh=finalsettings.r_thresh
  sievefac=finalsettings.sievefac
  cal_par=finalsettings.cal_par
  ;test a little
  
  ;Open input and ancillary files
  envi_open_file, infile, r_fid=fid,/no_realize,/no_interactive_query
  if (fid eq -1) then begin
    print,strtrim('Error opening input DWEL file',2)
    print,'Input File: '+strtrim(infile,2)
    err=1
    goto,out
  endif
  
  envi_open_file, ancfile, r_fid=ancillaryfile_fid,/no_realize,/no_interactive_query
  if (ancillaryfile_fid eq -1) then begin
    print,strtrim('Error opening input ancillary DWEL file',2)
    print,'Input File: '+strtrim(ancfile,2)
    err=2
    goto,out
  endif
  
  ; Get the file dimensions etc
  envi_file_query, fid, fname=infile, nb=nbands, nl=nlines, $
    ns=nsamples, bnames=bnames, wl=range, data_type=dt
  envi_file_query, ancillaryfile_fid, fname=ancfile, $
    nb=nb_anc, nl=nl_anc, ns=ns_anc, $
    data_type=anc_type, bnames=anc_bnames
    
  f_base=file_basename(infile)
  
  ;set up a base structure for the DWEL headers
  DWEL_headers={ $
    f_base:f_base $
    }
  ;find all of the DWEL headers in the hdr file as defined by FID
  status=dwel_get_headers(fid,DWEL_headers)
  
  ;Get date and time of the acquisition
  DWEL_date_time=''
  match = -1
  for i=0,n_elements(DWEL_headers.DWEL_scan_info)-1 do begin
    if (strmatch(DWEL_headers.DWEL_scan_info[i],'*DWEL Date Time*')) then match=i
  endfor
  if (match ge 0) then begin
    sf = strsplit(DWEL_headers.DWEL_scan_info[match],'=',/extract)
    DWEL_date_time = strtrim(strcompress(sf[1]),2)
  endif else begin
    DWEL_date_time = '2014-08-07-13-37'
  endelse
  DWEL_year=fix(strtrim(strmid(DWEL_date_time,0,4),2))
  
  print,'DWEL_year='+strtrim(string(DWEL_year),2)
  
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
  
  ;NOTE this will change! Both case in there - hope it persists
  ;Locate the DWEL Height
  DWEL_Height=-1.0
  match = -1
  for i=0,n_elements(DWEL_headers.DWEL_scan_info)-1 do begin
    if (strmatch(DWEL_headers.DWEL_scan_info[i],'*EVI Height*') or $
      strmatch(DWEL_headers.DWEL_scan_info[i],'*DWEL Height*')) then match=i
  endfor
  if (match ge 0) then begin
    sf = strsplit(DWEL_headers.DWEL_scan_info[match],'=',/extract)
    DWEL_Height = float(sf[1])
  endif else begin
    DWEL_Height = -1.0
  endelse
  
  print,'DWEL_Height='+strtrim(string(DWEL_Height),2)
  
  DWEL_pointcloud_info=[$
    'Program=DWEL_get_point_cloud',$
    'Title=Compute Point Cloud and Information Files',$
    'DWEL run Date & Time='+strtrim(DWEL_date_time,2),$
    'DWEL run Description='+strtrim(DWEL_description_record,2),$
    'DWEL Height='+strtrim(string(DWEL_height,format='(f10.3)'),2) $
    ]
    
  ; Check if file has been projected
  if (~DWEL_headers.proj_present) then begin
    projected = 0b
    projection_type='None'
    max_zen_ang=118.0
    mean_num_val=1
  endif else begin
    projected = 1b
    match = -1
    for i=0,n_elements(DWEL_headers.DWEL_projection_info)-1 do begin
      if (strmatch(DWEL_headers.DWEL_projection_info[i],'*Projection_type*')) then match=i
    endfor
    if (match ge 0) then begin
      sf = strsplit(DWEL_headers.DWEL_projection_info[match],'=',/extract)
      Projection_Type = strtrim(sf[1],2)
    endif else begin
      Projection_Type='None'
    endelse
  endelse
  
  if (projected) then begin
    match = -1
    for i=0,n_elements(DWEL_headers.DWEL_projection_info)-1 do begin
      if (strmatch(DWEL_headers.DWEL_projection_info[i],'*max_zenith_angle*')) then match=i
    endfor
    if (match ge 0) then begin
      sf = strsplit(DWEL_headers.DWEL_projection_info[match],'=',/extract)
      max_zen_ang = float(strtrim(sf[1],2))
    endif else begin
      max_zen_ang=118.0
    endelse

    match = -1
    for i=0,n_elements(DWEL_headers.DWEL_projection_info)-1 do begin
      if (strmatch(DWEL_headers.DWEL_projection_info[i],'*num_val_Stats*')) then match=i
    endfor
    if (match ge 0) then begin
      sf = strsplit(DWEL_headers.DWEL_projection_info[match],'=',/extract)
      sf = strmid(sf[1], 1, strlen(sf[1])-2)
      sf = strsplit(sf,',',/extract,count=tmpcount)
      mean_num_val = float(strtrim(sf[1], 2))
    endif else begin
      mean_num_val = 1
    endelse
  endif
  
  print,'max_zen_ang='+strtrim(string(max_zen_ang),2) 
  print, 'mean_num_val=', mean_num_val
 
  ; Check if file is apparent reflectance
  if (DWEL_headers.apprefl_present) then begin
    app_refl=1b
    cal_dat=0b
  endif else app_refl=0b
  
  ;Read the beam divergence from the scan info
  match = -1
  for i=0,n_elements(DWEL_headers.DWEL_scan_info)-1 do begin
    if (strmatch(DWEL_headers.DWEL_scan_info[i],'*Beam Divergence*')) then match=i
  endfor
  if (match ge 0) then begin
    sf = strtrim(strcompress(strsplit(DWEL_headers.DWEL_scan_info[match],'=',/extract)),2)
    t = strsplit(sf[1],' ',/extract)
    DWEL_div = float(t[0])
  endif else begin
    DWEL_div = -1.0
  endelse
  
  print,'Beam Divergence='+strtrim(string(DWEL_div),2)
  
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
    if (strpos(DWEL_headers.f_base,'1064') ne -1) then begin
      wavelength = 1064
    endif
    if (strpos(DWEL_headers.f_base,'1548') ne -1) then begin
      wavelength = 1548
    endif
  endelse
  
  print,'wavelength='+strtrim(string(wavelength),2)

  ; get the wire flag 
  base_info = DWEL_headers.DWEL_base_fix_info
  for i=0,n_elements(base_info)-1 do begin
    if (strmatch(base_info[i], '*Wire_Flag*', /fold_case)) then match=i
  endfor 
  if match ge 0 then begin
    text=strtrim(base_info[match],2)
    k=strpos(text,'=')
    wire_flag=fix(strtrim(strmid(text,k+1),2))
  endif else begin
    wire_flag = 0
  endelse 

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
    casing_type = 'LAM'
  endelse
  if (match ge 0) then print,'info match for casing type= ',strtrim(base_info[match],2)
  print,'casing type='+casing_type
  
  ;the calibration is now known as of 20141220
  tmp = where(cal_par le 0, tmpcount)
  if tmpcount gt 0 then begin
    ;; user does not provide calibration parameters or not provie a complete
    ;; correct one. use the default calibration paramters
    if (wavelength eq 1064) then begin
      if strcmp(laser_man, 'manlight', /fold_case) then begin
        if wire_flag then begin 
          ;; calibration parameters for scans with wire, but wire removed
          ;; NSF DWEL, for scaled intensity
          ;; estimated from stationary panel scans with wire
          ;; 10591.354 is the constant from unscaled return intensity DN
          ;; 2.305 is the nominal scale factor by onboard lambertian target
          ;; 0.984807753 is cosine(10 deg), panels were put 10 degrees slant
          ;; 1300 is the laser power setting of calibration scans.
          dwel_cal = 5863.906d0*0.984807753*5.727
          rpow = 1.402d0 ;; NSF DWEL
          eff_par = [3413.743d0, 0.895d0, 15.640d0, 3413.743d0]
        endif else begin
          ;; no-wire, same with wire-removed
          dwel_cal = 5863.906d0*0.984807753*5.727
          rpow = 1.402d0 ;; NSF DWEL
          eff_par = [3413.743d0, 0.895d0, 15.640d0, 3413.743d0]
        endelse 
      endif 
      if strcmp(laser_man, 'keopsys', /fold_case) then begin
        ;; from David, on 20141220
        dwel_cal = 2052936.432d0
        rpow = 1.9056d0
        eff_par = [6580.330d0, 0.3553d0, 43.396d0, 6580.330d0]
      endif 
    endif else begin
      if strcmp(laser_man, 'manlight', /fold_case) then begin
        if wire_flag then begin 
          ;; calibration parameters for scans with wire, but wire removed
          ;; NSF DWEL, for scaled intensity
          ;; estimated from stationary panel scans with wire
          ;; 10591.354 is the constant from unscaled return intensity DN
          ;; 2.305 is the nominal scale factor by onboard lambertian target
          ;; 0.984807753 is cosine(10 deg), panels were put 10 degrees slant
          ;; 1300 is the laser power setting of calibration scans.
          dwel_cal = 20543.960d0*0.984807753*5.103
          rpow = 1.566d0 ;; NSF DWEL
          eff_par = [5.133d0, 0.646d0, 1.114d0, 5.133d0]
        endif else begin
          ;; no-wire, same with wire-removed
          dwel_cal = 20543.960d0*0.984807753*5.103
          rpow = 1.566d0 ;; NSF DWEL
          eff_par = [5.133d0, 0.646d0, 1.114d0, 5.133d0]
        endelse 
      endif 
      if strcmp(laser_man, 'keopsys', /fold_case) then begin
        ;; from David, on 20141220
        dwel_cal = 712237.602d0
        rpow = 1.9056d0
        eff_par = [4483.089d0, 0.7317d0, 19.263d0, 4483.089d0]
      endif 
    endelse
  endif else begin
    ;; user has provided a valid calibration parameter set. 
    dwel_cal = cal_par[0]
    rpow = cal_par[5]
    eff_par = cal_par[1:4]
  endelse 
    
  ;set up the calibration as far as possible
  if (cal_dat) then begin
    s_Factor=1.0
    i_scale=1000.0
  endif else begin
    s_Factor=1.0
    i_scale=1.0
  endelse

  DWEL_pointcloud_info=[DWEL_pointcloud_info,$
    'DWEL beam wavelength='+strtrim(string(wavelength,format='(i10)'),2),$
    'DWEL beam divergence='+strtrim(string(DWEL_div,format='(f10.3)'),2),$
    'DWEL calibration Const='+strtrim(string(dwel_cal,format='(f10.3)'),2),$
    'DWEL calibration range power='+strtrim(string(rpow,format='(f10.3)'),2) $
    ]
  
  mask=bytarr(nsamples,nlines)
  zenith=fltarr(nsamples,nlines)
  azimuth=fltarr(nsamples,nlines)
  
  ; Get the Mask, Zenith, Azimuth from the ancillary file
  dims = [-1, 0, nsamples-1, 0, nlines-1]
  if (~projected) then begin
    Mask = byte(envi_get_data(fid=ancillaryfile_fid, dims=dims,pos=6))
    Zenith = float(envi_get_data(fid=ancillaryfile_fid, dims=dims,pos=7))/100.0
    Azimuth = float(envi_get_data(fid=ancillaryfile_fid, dims=dims,pos=8))/100.0
  endif else begin
    mask = byte(envi_get_data(fid=ancillaryfile_fid, dims=dims,pos=3))
    zenith = float(envi_get_data(fid=ancillaryfile_fid, dims=dims,pos=1))/100.0
    azimuth = float(envi_get_data(fid=ancillaryfile_fid, dims=dims,pos=2))/100.0
  endelse
  envi_file_mng,id=ancillaryfile_fid,/remove
  envi_file_mng, id=fid, /remove
  
  ; Call routine to set up pulse model
  dwel_itpulse_model_dual_nsf, wavelength, i_val, t_val, r_val, p_range, p_time, pulse, t_fwhm, r_fwhm
  
  ;S0 is a scale factor that relates the FWHM of the standard pulse and the mean FWHM of the data
  ;it has been chosen to be best for data in between ND015 (maybe 30% saturated) and ND100 (no saturated)
  ;
  S0=1.0  ;best for all but ND100
  fneg=-pulse[i_val[3]]
  h1=r_val[3]-r_val[2]
  h2=r_val[4]-r_val[2]
  
  print,'fneg='+strtrim(string(fneg),2)
  print,'h1='+strtrim(string(h1),2)
  print,'h2='+strtrim(string(h2),2)
  
  DWEL_pointcloud_info=[DWEL_pointcloud_info,$
    'S0='+strtrim(string(s0,format='(f12.4)'),2) $
    ]
    
  ; Resample pulse onto time base of data
  w = where((range ge p_range[i_val[0]]) and (range le p_range[i_val[6]]),numw)
  base_range=range[w[0]:w[numw-1]]/S0
  
  print,'number of elements in base range='+strtrim(string(n_elements(base_range)),2)
  print,'LH value='+strtrim(string(base_range[0]),2)
  print,'RH value='+strtrim(string(base_range[n_elements(base_range)-1]),2)
  
  p = interpol(pulse, p_range, base_range)
  
  value=min(abs(base_range-r_val[1]),min_ind)
  value=min(abs(base_range-r_val[5]),max_ind)
  
  range_step=base_range[1]-base_range[0]
  h1=round(h1/range_step)
  h2=round(h2/range_step)
  
  print,'range_step='+strtrim(string(range_step),2)
  print,'h1='+strtrim(string(h1),2)
  print,'h2='+strtrim(string(h2),2)
  
  ; Make sure we aren't asking for elements beyond the end of pulse array
  use_max_ind = max_ind < (n_elements(p)-1)
  if (use_max_ind ne max_ind) then begin
    print, 'Requested pulse array range truncated to element ',+strtrim(use_max_ind,2)
    print, ' '
  endif
  
  print,'min_ind='+strtrim(string(min_ind),2)
  ; Subset the pulse array
  max_ind=use_max_ind
  print,'max_ind='+strtrim(string(max_ind),2)
  p = p[min_ind:max_ind]
  
  ;options here are fixed for now
  if (app_refl) then begin
    threshold=0.005
    i_scale=1000.0
    DWEL_AppRefl=1b
  endif else begin
    threshold=9.5
    DWEL_AppRefl=0b
  endelse
  
  ;; now read in the right threshold from standard deviation of background noise
  ;; base! 
  fac=1.0
  if (DWEL_headers.filtfix_present) then begin
    match = -1
    for i=0,n_elements(DWEL_headers.dwel_filtered_fix_info)-1 do begin
      if (strmatch(DWEL_headers.dwel_filtered_fix_info[i],'*Noise_RMS=*')) then match=i
    endfor
    if (match ge 0) then begin
      sf = strsplit(DWEL_headers.dwel_filtered_fix_info[match],'=',/extract)
      threshold = float(strtrim(sf[1],2))
      print,'threshold from headers='+strtrim(string(threshold),2)
    endif
    match = -1
    for i=0,n_elements(DWEL_headers.dwel_filtered_fix_info)-1 do begin
      if (strmatch(DWEL_headers.dwel_filtered_fix_info[i],'*scale_mean=*')) then match=i
    endfor
    if (match ge 0) then begin
      sf = strsplit(DWEL_headers.dwel_filtered_fix_info[match],'=',/extract)
      scale_mean = float(strtrim(sf[1],2))
      print,'scale mean from headers='+strtrim(string(scale_mean),2)
    endif
  endif
  
  ;; read the average number of averaged shots in a projected bin

  ;; b_thresh=sdevfac*threshold
  ;; Because we have scaled waveforms in the filtered_fixbase processing
  ;; procedure, the standard deviation of background base level needs to be
  ;; scaled accordingly to reflect correct noise level and derive appropriate
  ;; threshold here.
  b_thresh=sdevfac*(scale_mean*threshold/sqrt(mean_num_val))

  sieve_thresh = sievefac*(scale_mean*threshold/sqrt(mean_num_val))
  
  print,'b_thresh='+strtrim(string(b_thresh),2)
  print, 'r_thresh=' + strtrim(string(r_thresh), 2)
  print, 'sieve_thresh=' + strtrim(string(sieve_thresh), 2)
  
  azimuth=azimuth-DWEL_az_n
  pos=where(azimuth lt 0.0,npos)
  if (npos gt 0) then azimuth[pos]=azimuth[pos]+360.0
  pos=0b
  pos=where(zenith gt max_zen_ang,num_zen)
  if (num_zen gt 0) then begin
    mask[pos]=0b
  endif
  pos=0b
  pos_mask=where(Mask ne 0b,num_mask)
  max_zenith=max(zenith[pos_mask])
  pos_mask=0b
  
  DWEL_pointcloud_info=[DWEL_pointcloud_info,$
    'DWEL_az_n='+strtrim(string(DWEL_az_n),2),$
    'B_Thresh='+strtrim(string(b_thresh),2), $
    'Sieve_Thresh='+strtrim(string(sieve_thresh),2) $
    ]
    
  ; ***********************
  ; Get path and file name as separate strings
    
  last=strpos(infile,path_sep(),/reverse_search)
  in_path=file_dirname(infile)
  in_base=strtrim(strmid(infile,last+1,strlen(infile)-last-1),2)
  
  last=strpos(ancfile,path_sep(),/reverse_search)
  anc_path = file_dirname(ancfile)
  anc_base=strtrim(strmid(ancfile,last+1,strlen(ancfile)-last-1),2)
  
  ;Set up the metadata elements
  
  Processing_Date_Time=''
  Run_Number=runcode
  Description=DWEL_description_record
  Input_Path=in_path
  Input_File=in_base
  Acquisition_Date_Time=DWEL_Date_Time
  Ancillary_Path=anc_path
  Ancillary_File=anc_base
  Projection=Projection_Type
  DWEL_AppRefl=DWEL_AppRefl
  DWEL_Height=DWEL_Height
  DWEL_Az_North=DWEL_az_n
  Max_Zenith_Angle=Max_Zenith
  Range_Step=abs(range[1]-range[0])
  Threshold=threshold
  b_thresh=b_thresh
  Zero_Hit_Option=save_zero_hits
  Add_DWEL=Add_DWEL
  X_scale=1.0
  Y_scale=1.0
  Z_scale=1.0
  X_offset=0.0
  Y_offset=0.0
  Z_offset=0.0
  I_Scale=i_scale
  Point_File_Path=''
  Point_File_name=''
  Zero_Hit_Number=0L
  Shot_Hit_Number=0L
  Total_Hit_Number=0L
  Nrecs=0L
  Max_X=0.0
  Min_X=0.0
  Max_Y=0.0
  Min_Y=0.0
  Max_Z=0.0
  Min_Z=0.0
  Min_Intensity=0.0
  Max_Intensity=0.0
  
  ;***************
  ;now set up structures and push onto the stack
  
  meta_data={pbmeta, $
    Processing_Date_Time:Processing_Date_Time,$
    Run_Number:Run_Number,$
    Description:Description,$
    Input_Path:Input_Path,$
    Input_File:Input_File,$
    Acquisition_Date_Time:Acquisition_Date_Time,$
    Ancillary_Path:Ancillary_Path,$
    Ancillary_File:Ancillary_File,$
    Projection:Projection,$
    DWEL_AppRefl:DWEL_AppRefl,$
    DWEL_Height:DWEL_Height,$
    DWEL_Az_North:DWEL_Az_North,$
    Max_Zenith_Angle:Max_Zenith_Angle,$
    Range_Step:Range_Step,$
    Threshold:Threshold,$
    b_thresh:b_thresh,$
    r_thresh:r_thresh,$
    sieve_thresh:sieve_thresh,$
    Fneg:fneg,$
    h1:h1,$
    h2:h2,$
    xmin:xmin,$
    xmax:xmax,$
    ymin:ymin,$
    ymax:ymax,$
    zlow:zlow,$
    zhigh:zhigh,$
    wavelength:wavelength,$
    Zero_Hit_Option:Zero_Hit_Option,$
    Add_DWEL:Add_DWEL,$
    save_br:save_br,$
    save_pfilt:save_pfilt,$
    X_scale:X_scale,$
    Y_scale:Y_scale,$
    Z_scale:Z_scale,$
    X_offset:X_offset,$
    Y_offset:Y_offset,$
    Z_offset:Z_offset,$
    I_Scale:I_Scale,$
    Point_File_Path:point_file_Path,$
    Point_File_name:point_file_name,$
    Zero_Hit_Number:Zero_Hit_Number,$
    Shot_Hit_Number:Shot_Hit_Number,$
    Total_Hit_Number:Total_Hit_Number,$
    Nrecs:Nrecs,$
    Max_X:Max_X,$
    Min_X:Min_X,$
    Max_Y:Max_Y,$
    Min_Y:Min_Y,$
    Max_Z:Max_Z,$
    Min_Z:Min_Z,$
    Min_Intensity:Min_Intensity,$
    Max_Intensity:Max_Intensity,$
    laser_man:laser_man,$
    wire_flag:wire_flag $
    }
    
  pb_meta=ptr_new(meta_data,/no_copy)
  
  out_file=strtrim(outfile,2)
  n_base=strlen(out_file)
  n_dot=strpos(out_file,'.',/reverse_search)
  if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
    oi_name=out_file+'_pcinfo.img'
  endif else begin
    oi_name=strmid(out_file,0,n_dot)+'_pcinfo.img'
  endelse
  oi_name=strtrim(oi_name,2)
  
  ;see if the output image file exists & remove if it does!
  if(file_test(oi_name)) then begin
    fids=envi_get_file_ids()
    if(fids[0] eq -1) then begin
      file_delete, oi_name,/quiet
      print,'old image file deleted'
    endif else begin
      for i=0,n_elements(fids)-1 do begin
        envi_file_query,fids[i],fname=tname
        if (strtrim(strlowcase(oi_name),2) eq $
          strtrim(strlowcase(tname),2)) then begin
          envi_file_mng,id=fids[i],/remove
          print,'old image file removed from ENVI'
        endif
      endfor
      file_delete, oi_name,/quiet
      print,'old image file deleted'
    endelse
  endif
  
  base_stats={pbstats, $
    infile:infile,$
    outfile:out_file,$
    oi_name:oi_name,$
    range:range,$
    mask:mask,$
    zenith:zenith,$
    azimuth:azimuth,$
    cal_dat:cal_dat,$
    DWEL_div:DWEL_div,$
    Dwel_cal:Dwel_cal,$
    rpow:rpow,$
    eff_par:eff_par,$
    s_Factor:s_Factor,$
    pulse:pulse,$
    p_range:p_range}
    
  pb_stats=ptr_new(base_stats,/no_copy)
  pb_info = ''
  
  ;Apply the point cloud processing
  ;********************************************
  perr=0b
  dwel_apply_ptcl_filter,p,pb_stats,pb_meta,pb_info,error=perr
  ;********************************************
  
  if (perr ne 0b) then begin
    print,'Error in point cloud processing'
    print,'Error called from apply_ptcl_filter'
    print,'Local error number='+strtrim(string(perr),2)
    ptr_free,pb_stats
    ptr_free,pb_meta
    err=99
    goto,out
  endif
  
  DWEL_pointcloud_info=[DWEL_pointcloud_info,pb_info]
  pb_info=''
  DWEL_pointcloud_info=[DWEL_pointcloud_info,$
    'PointCloud Info File='+strtrim(oi_name,2),$
    'Output PointCloud File='+strtrim(out_file,2) $
    ]
    
  tname=(*pb_stats).oi_name
  ;write out header for info image
  descrip='Info image for Point Cloud '+strtrim(out_file,2)
  bnames=['Nhits','Total d','Total d0','Total I','Total I2','First_Hit','Mean_Range*100','Last_Hit','Zenith*100','Azimuth*100','Residual','Gap','Mask','Mean_z']
  envi_setup_head,fname=tname,ns=nsamples,nl=nlines,nb=14,$
    xstart=0,ystart=0,$
    data_type=3, interleave=0, bnames=bnames, $
    descrip=descrip, /write
    
  envi_open_file,(*pb_stats).oi_name,r_fid=anc_fid,/no_interactive_query,/no_realize
  ;write out the previous header records
  status=dwel_put_headers(anc_fid,DWEL_headers)
  envi_assign_header_value, fid=anc_fid, keyword='DWEL_pointcloud_info', $
    value=DWEL_pointcloud_info
    
  envi_write_file_header, anc_fid
  envi_file_mng,id=anc_fid,/remove
  dwel_pointcloud_info=['']
  
  ptr_free,pb_stats
  ptr_free,pb_meta
  
  ; **************
  ; Remove all remaining file handles from ENVI
  envi_file_mng, id=ancillaryfile_fid, /remove
  envi_file_mng, id=fid, /remove
  ancillaryfile_fid=0b
  fid=0b
  if(ptr_valid(pb_stats)) then ptr_free,pb_stats
  if(ptr_valid(pb_meta)) then ptr_free,pb_meta
  if(ptr_valid(p_stat)) then ptr_free,p_stat
  pb_stats=0b
  pb_meta=0b
  p_stat=0b
  ; **************
  print,'Completed point_cloud'
  print,'Output PointCloud File Base='+strtrim(out_file,2)
  print,'PointCloud Info File='+strtrim(oi_name,2)
  
  print,' '
  print,'*************************************'
  print,' '
  
  out:
  if (ptr_valid(pb_stats)) then ptr_free,pb_stats
  if (ptr_valid(pb_meta)) then ptr_free,pb_meta
  if(ptr_valid(p_stat)) then ptr_free,p_stat
  free_lun,tfile,/force
  
  heap_gc,/verbose

  ;; write processing time summary
  print, '*****************************************'
  print, 'Processing program = dwel_get_point_cloud'
  print, 'Input cube file size = ' + $
    strtrim(string(double(procfilesize)/(1024.0*1024.0*1024.0)), 2) + ' G'
  print, 'Processing time = ' + strtrim(string((systime(1) - starttime)), $
    2) + ' ' + $
    'seconds'
  print, '*****************************************'
  
  return
  
end
;======================================================================
