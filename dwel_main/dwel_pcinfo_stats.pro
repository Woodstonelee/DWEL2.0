;======================================================================
;  The following set of procedures is part of the CSIRO EOC Hyperion
;  Workshop module set. The routines all have a similar structure.
;  They can be used stand-alone or as part of a Project structure with
;  ENVI main menu button.
;
;  The general structure for a procedure called "name" is:
;
;  [Procedures called from name-doit]
;    These may be called anything
;    they arerranged so that every procedure occurs before its first
;    reference. This is good IDL practice.
;  pro name_doit
;    the main area of action for the algorithm - the DOIT!
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
pro DWEL_pcinfo_stats_doit, event
  compile_opt idl2
  
  outfile=''
  odfile=77
  
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
    result=dialog_message(info_text,/error,title='Error in DWEL_pcinfo_stats')
    ;now see if the output file has a handle
    goto, cleanup
  endif
  
  ;get the grandfather widget in case...
  widget_control,event.top,get_uvalue=pstate
  
  ;clean up any fids which are no longer where they were!
  ;ENVI issue that is annoying and leads to confusion
  clean_envi_file_fids
  
  ;define the input and output images as strings
  image=''
  cal_name=''
  out_name=''
  o_name=''
  ofile=99
  
  d_name=''
  out_dname=''
  od_name=''
  odfile=77
  
  ;select the input image file
  ;the band selection is turned off as it is done in the ascii file
  input:
  envi_select, dims=dims,fid=fid,title='Select DWEL pc_info image for Stats', $
    pos=pos_b, /file_only, /no_spec
  if(fid eq -1) then begin
    result=dialog_message('Try again ? (No/Yes) ',/question,/default_no, $
      title='No file selected to fix or operation cancelled !',$
      dialog_parent=event.top)
    if(result eq 'No') then begin
      goto, cleanup
    endif else goto, input
  endif
  
  ;get the input image dimensions and other info
  envi_file_query, fid, ns=ns, nl=nl, nb=nb, fname=image, $
    byte_swap=order, data_type=type, interleave=ftype, $
    bnames=band_name, $
    xstart=xstart, ystart=ystart,wl=wl
    
  samples=ns
  lines=nl
  bands=nb
  
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
  last=strpos(image,path_sep(),/reverse_search)
  f_path=strmid(image,0,last+1)
  f_base=strtrim(strmid(image,last+1,strlen(image)-last-1),2)
  
  n_base=strlen(image)
  n_dot=strpos(image,'.',/reverse_search)
  
  ;Set up the List the image information to check
  Info_Text=[$
    'Base information for: '+strtrim(f_base,2),$
    'Path is: '+strtrim(f_path,2),$
    'File interleave is: '+strtrim(ft_nam,2),$
    'Data type is: '+strtrim(dt_nam,2),$
    ' ',$
    'samples in image: '+strtrim(string(samples,format='(i5)'),2),$
    'lines in image: '+strtrim(string(lines,format='(i5)'),2),$
    'bands in image: '+strtrim(string(bands,format='(i5)'),2),$
    ' ',$
    ' ',$
    'All OK ? (Y/N)'$
    ]
    
  ;now get the DWEL headers that are present
  ;set up a base structure for the DWEL headers
  DWEL_headers={ $
    f_base:f_base $
    }
    
  ;find all of the DWEL headers in the hdr file as defined by FID
  status=dwel_get_headers(fid,DWEL_headers)
  
  if (not status) then begin
    result=dialog_message('Bad FID in get_headers - Restart ? ',$
      /question, title='DWEL Header setup cancelled', $
      dialog_parent=event.top)
    if(result eq 'No') then begin
      goto, cleanup
    endif else goto, input
  endif
  
  if (DWEL_headers.headers_present le 0s or ~DWEL_headers.run_present) then begin
    result=dialog_message('NOT a valid DWEL file - Restart ? ',$
      /question, title='Invalid DWEL file attached', $
      dialog_parent=event.top)
    if(result eq 'No') then begin
      goto, cleanup
    endif else goto, input
  endif
  
  if (~DWEL_headers.ptcld_present) then begin
    print,'bad file, not ptcld header info'
    info_note=[ $
      'File header indicates DWEL data NOT from point cloud processing.',$
      '',$
      'Try another input file? (No/Yes)' $
      ]
    result=dialog_message(info_note,/question,$
      title='Bad input file')
    if(result eq 'Yes') then goto, input else goto, cleanup
  endif
  
  ;Get the Total_d stats
  buf=''
  match = -1
  for i=0,n_elements(DWEL_headers.DWEL_pointcloud_info)-1 do begin
    if (strmatch(DWEL_headers.DWEL_pointcloud_info[i],'*d_Stats*') $
      or strmatch(DWEL_headers.DWEL_pointcloud_info[i],'*d_Stats*')) then match=i
  endfor
  if (match ge 0) then begin
    print,'match=',DWEL_headers.DWEL_pointcloud_info[match]
    sf = strsplit(DWEL_headers.DWEL_pointcloud_info[match],'=',/extract)
    if (n_elements(sf) gt 1) then begin
      buf = strtrim(strcompress(sf[1]),2)
    endif else begin
      buf = ''
    endelse
  endif else begin
    buf = ''
  endelse
  if ((strpos(buf,'(') eq 0) and (strpos(buf,')',/reverse_search) eq strlen(buf)-1)) then begin
    buf=strmid(buf,1,strlen(buf)-2)
  endif
  val=float(strsplit(buf,',',/extract))
  
  d_min=val[0]
  d_max=val[2]
  
  for j=0,n_elements(DWEL_headers.DWEL_pointcloud_info)-1 do begin
    print,DWEL_headers.DWEL_pointcloud_info[j]
  endfor
  
  ;set info into prompt and display it
  Info_Text=[Info_Text[0:8],DWEL_headers.DWEL_pointcloud_info,Info_Text[8:10]]
  ;if all is OK get on with it!
  result=dialog_message(Info_Text,/question, $
    title='OK to get on with DWEL Pgap Estimation ? (Y/N)', $
    dialog_parent=event.top)
  if (result eq 'No') then goto, input
  
  ;-------------------------------------------------------------------------
  
  ;now get the output Stats file name
  
  output_dist:
  
  n_base=strlen(f_base)
  n_dot=strpos(f_base,'.',/reverse_search)
  
  if (d_name eq '') then begin
    if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
      od_name=f_base+'_pc_info_stats.txt'
    endif else begin
      od_name=strmid(f_base,0,n_dot)+'_pc_info_stats.txt'
    endelse
  endif
  
  out_dname=dialog_pickfile(title='Select pc_info Stats file name',file=od_name, $
    dialog_parent=event.top,path=f_path)
    
  ;now check for no name, output=input or existing file
  if(out_dname eq '') then begin
    result=dialog_message('Try again ? (No/Yes) ',/question,/default_no, $
      title='No pc_info Stats file selected!',$
      dialog_parent=event.top)
    if(result eq 'No') then begin
      goto, cleanup
    endif else goto, output_dist
  ;We will not let the output be the input file
  endif else if (strtrim(out_dname,2) eq $
    strtrim(image,2)) then begin
    info_text=[strtrim('pc_info Stats File name conflict',2),$
      strtrim('Output cannot be the Input !',2),$
      strtrim('Try another selection',2)]
    result=dialog_message(info_text,/error,title='Invalid Input',$
      dialog_parent=event.top)
    out_dname=''
    od_name=''
    goto, output_dist
  endif
  
  ;see if the output file exists
  if(file_test(out_dname)) then begin
    result=dialog_message('OK to overwrite '+strtrim(out_dname,2)+' ? (Yes/No) ', $
      /question, $
      title='Output pc_info Stats File Exists',$
      dialog_parent=event.top)
    if(result eq 'Yes') then begin
      fids=envi_get_file_ids()
      if(fids[0] eq -1) then begin
        file_delete,out_dname,/quiet
      endif else begin
        for i=0,n_elements(fids)-1 do begin
          envi_file_query,fids[i],fname=tname
          if (strtrim(strlowcase(out_dname),2) eq $
            strtrim(strlowcase(tname),2)) then begin
            envi_file_mng,id=fids[i],/remove
          endif
        endfor
        file_delete,out_dname,/quiet
      endelse
    endif else begin
      out_dname=''
      od_name=''
      goto, output_dist
    endelse
  endif
  ;open the output file
  openw,odfile,out_dname,/get_lun,error=err
  if (err ne 0) then begin
    print,'Error opening the pcinfo output file'
    print,'File name: '+strtrim(out_dname,2)
    print,!error_state.msg
    goto,cleanup
  endif
  printf,odfile,strtrim('Process PCInfo File to get G',2)
  printf,odfile,strtrim('Run made at: '+systime(),2)
  printf,odfile,'Input File='+strtrim(image,2)
  flush,odfile
  
  ;-------------------------------------------------------------------------
  
  compute_Stats:
  
  dfac=(d_max-d_min)/4095.0
  
  Mask=envi_get_data(fid=fid,dims=dims,pos=10)
  Nhits=envi_get_data(fid=fid,dims=dims,pos=0)
  
  printf,odfile,'Stats Format=Min,Mean,Max,Sdev'
  ;get nhits stats
  image_statistics,Nhits,mask=mask,minimum=emin,maximum=emax,mean=emean,stddev=esdev
  printf,odfile,'Nhits Stats='+strtrim(strcompress(strjoin(string([emin,emean,emax,esdev]),','),/remove_all),2)
  
  total_d=dfac*float(envi_get_data(fid=fid,dims=dims,pos=1))+d_min
  
  ;get total_d stats
  image_statistics,total_d,mask=mask,minimum=emin,maximum=emax,mean=emean,stddev=esdev
  printf,odfile,'Total_d Stats='+strtrim(strcompress(strjoin(string([emin,emean,emax,esdev]),','),/remove_all),2)
  
  numbins=40
  hmax=1.5
  step=hmax/float(numbins)
  
  pos=where((mask ne 0)and (Nhits eq 1),npos)
  if (npos gt 0) then $
    hist=histogram(total_d[pos],min=0.0,max=hmax,nbins=numbins) $
  else $
    hist=histogram(total_d,min=0.0,max=hmax,nbins=numbins)
    
  bin_ind=findgen(numbins)*step
  bufx=strtrim(string(bin_ind),2)
  bufx=strtrim(strcompress(strjoin(bufx,','),/remove_all),2)
  
  printf,odfile,'d histogram'
  printf,odfile,bufx
  buf=strtrim(string(hist),2)
  buf=strtrim(strcompress(strjoin(buf,','),/remove_all),2)
  printf,odfile,buf
  flush,odfile
  
  range=float(envi_get_data(fid=fid,dims=dims,pos=5))/100.0
  ;get range stats
  image_statistics,range,mask=mask,minimum=emin,maximum=emax,mean=emean,stddev=esdev
  print,'range Stats=',emin,emean,emax,esdev
  
  zenith=float(envi_get_data(fid=fid,dims=dims,pos=6))/10.0
  
  step=5.0
  nstep=24
  
  printf,odfile,'j,Zenith,Number,Min,Mean,Max,SDev,99.9%'
  flush,odfile
  
  v99=0.999
  zth=fltarr(nstep)
  yth=fltarr(nstep)
  
  num_reg=0
  
  for j=0,nstep-1 do begin
    ;
    zmin=0.0+float(j)*step
    zmax=zmin+step
    pos=where((mask ne 0)and (Nhits eq 1) and (zenith gt zmin) and (zenith le zmax),numpos)
    if (numpos le 0) then begin
      printf,odfile,strtrim(strcompress(strjoin(string([float(j+1),(zmin+zmax)/2.0,float(numpos),0.0,0.0,0.0,0.0,0.0]),','),/remove_all),2)
      goto,nzenith
    endif
    tmask=intarr(ns,nl)
    if (numpos gt 0) then tmask[pos]=1
    image_statistics,total_d,mask=tmask,minimum=emin,maximum=emax,mean=emean,stddev=esdev
    hmin=emin
    hmax=emax
    hist=histogram(total_d[pos],min=hmin,max=hmax,nbins=numbins)
    
    div=total(hist)
    if (div ne numpos) then print,'WARNING! div ne numpos',j,div,numpos
    ;
    cumulative=total(float(hist),/cum)/float(div)
    posb=where(cumulative gt v99,numposb)
    if (numposb le 0) then begin
      print,'numposb le 0! at ',j
    endif
    ;
    ind1=posb[0]-1
    ind2=ind1+1
    
    if (ind1 lt 0) then begin
      est99=emax
    endif else begin
      if ((cumulative[ind1] gt v99) or (cumulative[ind2] lt v99)) then begin
        print,'bad bracket at j=',j
        print,'ind1=',ind1,' ind2=',ind2
      endif
      est99=(bin_ind[ind2]*(v99-cumulative[ind1])+bin_ind[ind1]*(cumulative[ind2]-v99))/ $
        (cumulative[ind2]-cumulative[ind1])
    endelse
    
    printf,odfile,strtrim(strcompress(strjoin(string([float(j+1),(zmin+zmax)/2.0,float(numpos),emin,emean,emax,esdev,est99]),','),/remove_all),2)
    ;
    flush,odfile
    zth[num_reg]=(zmin+zmax)/2.0
    yth[num_reg]=est99
    num_reg=num_reg+1
    pos=0b
    tmask=0b
    
    nzenith:
  endfor
  
  help,zth
  help,yth
  t1=cos(!dtor*zth)
  t2=(2/!pi)*sin(!dtor*zth)
  x1=t1##transpose(t1)
  x2=t1##transpose(t2)
  x3=t2##transpose(t2)
  
  h1=t1##transpose(yth)
  h2=t2##transpose(yth)
  
  delta=x1*x3-x2^2
  
  if (abs(delta) lt 0.00001) then begin
    print,'delta='+strtrim(string(delta),2)
    print,'delta is TOO small - expect problems!'
  endif
  ah=(x3*h1-x2*h2)/delta
  av=(-x2*h1+x1*h2)/delta
  print,'ah,av='+strtrim(string(ah),2)+','+strtrim(string(av),2)
  a_star=ah+av
  fv=av/a_star
  
  if ((fv lt 0.0) or (fv gt 1.0)) then begin
    print,'fv=',fv
    print,'fv is NEGATIVE or gt 1 - problems!'
  endif
  
  model=ah*t1+av*t2
  err=sqrt(total((yth-model)^2)/float(num_reg))
  
  print,'err='+strtrim(string(err),2)
  
  printf,odfile,''
  printf,odfile,'Model Output'
  printf,odfile,'Model Num='+strtrim(string(num_reg),2)
  printf,odfile,'Model A='+strtrim(string(a_star),2)
  printf,odfile,'Model Fv='+strtrim(string(fv),2)
  printf,odfile,'Model RMS='+strtrim(string(err),2)
  
  flush,odfile
  free_lun, odfile,/force
  
  print,'dwel_pcinfo_stats finished successfully'
  
  ;exit to caller
  cleanup:
  pout=0b
  band_pos=0b
  gain=0b
  offset=0b
  
  print,'exiting dwel_pcinfo_stats'
  
  free_lun, odfile,/force
  ;get state pointer
  widget_control,event.top,get_uvalue=pstate
  ;clean up pointers
  widget_control,event.top,/destroy
  
  heap_gc,/verbose
  return
  
end

;======================================================================

function DWEL_pcinfo_stats_bannertext
  compile_opt idl2
  
  ;Enter the banner text as a string array
  
  text=[ $
    'Module Description:',$
    ' ',$
    'Procedure DWEL_pcinfo_stats computes a set of statistics to provide',$
    'estimates of pure gaps and hard and soft targets in DWEL data',$
    'The statistics are mainly based on treating the returns as a',$
    'distribution or histogram in range. Pure gaps have effectively no',$
    'information, hard targets have single hard returns with the shape of',$
    'the outgoing pulse and soft targets have distributed returns.',$
    'The statistics are normally carried out on Apparent Reflectance',$
    'processed images',$
    'The data can be selected and combined with other statistics',$
    'to see if soft and hard targets and other components of the forest',$
    'can be classified or recognised' $
    ]
    
  return,text
  
end

pro DWEL_pcinfo_stats_action_go, event
  compile_opt idl2
  
  ;Event handler for hitting the "GO" button
  ;simply calls the main action routine
  
  widget_control,event.top,get_uvalue=pstate
  
  (*pstate).go=1
  
  DWEL_pcinfo_stats_doit, event
  
  return
  
end

pro DWEL_pcinfo_stats_action_exit, event
  compile_opt idl2
  
  ;Event handler for hitting the "EXIT" button
  ;cleans up and destroys the widget hierarchy
  
  widget_control,event.top,get_uvalue=pstate
  ptr_free, pstate
  widget_control,event.top,/destroy
  
  ;some protective programming for now
  heap_gc,/verbose
  
  return
  
end

pro DWEL_pcinfo_stats_resize, event
  compile_opt idl2
  
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
    result=dialog_message(info_text,/error,title='Error in DWEL_pcinfo_stats_resize', $
      dialog_parent=event.top)
    goto, cleanup
  endif
  
  ;Window resize event handler
  
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

pro DWEL_pcinfo_stats_cleanup,top
  compile_opt idl2
  
  ;The cleanup routine for Xmanager
  
  ;get state pointer
  widget_control,top,get_uvalue=pstate
  
  ;clean up pointers
  ptr_free, pstate
  
  return
  
end

pro DWEL_pcinfo_stats, event
  compile_opt idl2
  
  ;The main routine that sets up the widget hierarchy and manages the events
  ;There are only two main events here so it is not so hard
  
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
    result=dialog_message(info_text,/error,title='Error in DWEL_pcinfo_stats')
    goto, cleanup
  endif
  
  xoff=200
  yoff=150
  
  ;now set up the top level base widget
  envi_center, xoff,yoff
  wb_base=widget_base(title='Compute Gap Stats',/column,$
    xoffset=xoff,yoffset=yoff,/frame,$
    /tlb_size_events)
    
  rev=1s
  
  ;get the procedure description and disclaimer to display
  intro_text=DWEL_pcinfo_stats_bannertext()
  dis_claim=DWEL_disclaimer(rev)
  
  ;Set up a base widget to hold the description and disclaimer
  wb_info=widget_base(wb_base,/column)
  
  ;Use two slider label box routines from the ENVI library
  w_2=widget_slabel(wb_info,prompt=intro_text,/frame,xsize=53,ysize=10)
  w_3=widget_slabel(wb_info,prompt=dis_claim,/frame,xsize=53,ysize=1)
  
  ;Now set up another base widget for the action buttons
  wb_action=widget_base(wb_base,/row)
  
  ;Action buttons are just "go" and "exit"
  w_go=widget_button(wb_action,value='Go',uvalue='DWEL_pcinfo_stats_action_go',$
    event_pro='DWEL_pcinfo_stats_action_go',/dynamic_resize)
  w_exit=widget_button(wb_action,Value='Exit',uvalue='DWEL_pcinfo_stats_action_exit',$
    event_pro='DWEL_pcinfo_stats_action_exit',/dynamic_resize)
    
  ;realise the widget hierarchy
  widget_control,wb_base,/realize
  
  go=0
  
  widget_geometry_1=widget_info(wb_base,/geometry)
  widget_geometry_2=widget_info(wb_action,/geometry)
  
  xy_old=[widget_geometry_1.xsize,widget_geometry_1.ysize]
  xy_act=[widget_geometry_2.xsize,widget_geometry_2.ysize]
  
  ;set up the state structure
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
  xmanager,'DWEL_pcinfo_stats',wb_base, $
    event_handler='DWEL_pcinfo_stats_resize',$
    cleanup='DWEL_pcinfo_stats_cleanup',/no_block
    
  ;exit back to ENVI
  ;leave a bit of a trace for bad exits for now
  cleanup:
  
  heap_gc,/verbose
  
  return
  
end
