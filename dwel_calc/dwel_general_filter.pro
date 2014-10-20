;
;+
; NAME:
;DWEL_GENERAL_FILTER
;
; PURPOSE:
;Convolve the pulse model with waveforms. In such way each return pulse will
;become an iterated pulse with symmetrical side lobes and some noise may be
;filtered. But the transient ringing noise could be enhanced. 
;
; CATEGORY:
;DWEL waveform processing, pulse filtering.
;
; CALLING SEQUENCE:
;dwel_general_filter, infile, p, mpos, outfile, ierr
;
; INPUTS:
;infile = input file name of DWEL cube to be filtered.
;p = the pulse model to filter waveforms. 
;mpos = the origin (max value) point of the filter, usually the centre, the
;index to this location.
;outfile = output file name of the filtered DWEL waveform cube. 
; ierr = a variable to receive error code.
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
;Convolve the pulse model with each waveform scan line by scan line and write
;the convolution results to the output file.  
;
; MODIFICATION HISTORY:
;David Jupp, Sept 2014 - Created this routine. 
;Zhan Li, Oct 2014 - Added documentation comments.
;-
function dwel_general_filter, infile, p, mpos, outfile, ierr
  compile_opt idl2
  envi, /restore_base_save_files
  envi_batch_init, /no_status_window
  ;
  ;general 1D filter on bands program
  ;NOTE: assumes file is BIL, if not it returns error
  ierr=0
  ofile=30
  inlun=35
  if ((n_elements(p) le 0) or ((mpos lt 0) $
    or (mpos gt n_elements(p)-1))) then begin
    ierr=1
    goto,cleanup
  endif
  ;here open the input file and get info eg:
  ;  fid, dt, nsamples, nlines, nbands
  ;
  envi_open_file, infile, r_fid=fid, /no_interactive_query, /no_realize
  if(fid eq -1) then begin
    ierr=2
    goto, cleanup
  endif
  ;get the input image dimensions and other info
  envi_file_query, fid, ns=ns, nl=nl, nb=nb, $
    data_type=dt, interleave=ftype, wl=wl,$
    dims=dims, file_type=f_type
    
  result=envi_file_type(f_type)
  ;  print,'input file type = ',result
  if (result ne 'ENVI Standard') then begin
    ierr=3
    envi_file_mng,id=fid,/remove
    goto, cleanup
  endif
  ;
  if (ftype ne 1) then begin
    ierr=4
    envi_file_mng,id=fid,/remove
    goto, cleanup
  endif
  
  nsamples=ns
  nlines=nl
  nbands=nb
  
  nbytes=dt2nb(dt) ; get the number of bytes from the data type. 
  
  ;set data type (based on input file)
  if (dt lt 4 or dt gt 9) then begin
    zero=0s
  endif else zero=0.0
  
  ;close up the envi file
  envi_file_mng,id=fid,/remove
  
  ; Calculate some values needed in the convolution
  ; mpos is the origin (max value) point of the filter - usually the centre
  
  nlead = mpos
  ntrail = n_elements(p)-nlead-1
  
  ; Pad pulse with zeros to make symmetrical about the peak
  ppad=float(p)
  if (ntrail gt nlead) then begin
    ppad = [replicate(0.0,ntrail-nlead),p]
  endif else if (nlead gt ntrail) then begin
    ppad = [p,replicate(0.0,nlead-ntrail)]
  endif
  value=max(ppad,max_point)
  
  ;print,'elements in ppad=',n_elements(ppad)
  ;print,'mid point=',n_elements(ppad)/2
  print,'max_point=',max_point
  print,'max_value=',value
  
  num_pad=n_elements(ppad)
  fnum=float(num_pad)
  psum=total(ppad)
  pp_filter=fltarr(1,num_pad)
  pp_filter[0,*]=float(ppad)
  ppad=0b
  
  ;  print,'num_pad=',num_pad
  ;  print,'psum=',psum
  
  err=0
  ;open input file
  openr,inlun,infile,/get_lun,error=err
  if (err gt 0) then begin
    ierr=5
    goto,cleanup
  endif
  err=0
  ; Open output file
  openw,ofile,outfile,/get_lun,error=err
  if (err ne 0) then begin
    ierr=6
    goto, cleanup
  endif
  
  bufrs=long64(nbytes)*long64(nbands)*long64(nsamples)
  pointsz=long64(0)
  
  pos=indgen(nbands)
  ; Loop through each line, apply filter.
  ; then Write to output file.
  for i=0L,nlines-1L do begin
    ;    line = envi_get_slice(fid=fid, /bil, pos=pos, line=i, xs=0, xe=nsamples-1)  ; line is ns by nb
    line=read_binary(inlun, data_start=pointsz, data_dims=[nsamples,nbands], $
      data_type=dt)
    pointsz=long64(pointsz) + long64(bufrs)
    b=fltarr(nsamples,nbands)
    ;
    T_Pad=total(float(line[*,0:29]),2)/30.0 ; ns (column) by 1 (row) array, each
                                ; element is an average of the first 30 bins of
                                ; each waveform
    T_pad=cmreplicate(T_pad,num_pad)        ; ns (column) by num_pad (row)
                                ; array, each row is the averages of all columns
                                ; in this line 
    ;    if (i eq 0) then help,t_pad
    B_pad=total(float(line[*,nbands-30:nbands-1]),2)/30.0 ; average of the last 30 bins
    B_pad=cmreplicate(B_pad,num_pad)
    ;    if (i eq 0) then help,b_pad
    t=[[T_pad],[float(line)],[B_pad]] ; use the average of first and last 30 bins to pad the leading and trail part of all the waveforms in this line
    ;    if (i eq 0) then help,t
    T_pad=0b
    B_pad=0b
    line=0b
    c=convol(t,pp_filter,psum)
    ;    if (i eq 0) then help,c
    b=reform(c[*,num_pad:num_pad+nbands-1])
    ;    if (i eq 0) then help,b
    c=0b
    t=0b
    ;
    writeu, ofile, fix(round(b))
    b=0b
  endfor
  free_lun, ofile,/force
  ;
  dt_br=2 ; integer
  envi_open_file, infile,r_fid=fid,/no_interactive_query,/no_realize
  
  ;get the input image dimensions and other info
  envi_file_query, fid, ns=ns, nl=nl, nb=nb, $
    data_type=dt, interleave=ftype, wl=wl,$
    dims=dims, file_type=f_type
  inherit=envi_set_inheritance(fid,dims,pos,/full)
  
  descrip='DWEL_Oz general filter applied to '+strtrim(infile,2)
  envi_setup_head,fname=outfile,ns=nsamples,nl=nlines,nb=nbands,$
    xstart=0,ystart=0,inherit=inherit,wl=wl,$
    data_type=dt_br, interleave=1, descrip=descrip, /write
    
  envi_file_mng,id=fid,/remove
  
  cleanup:
  ; Close output files in case
  free_lun, ofile,/force
  ;
  return, ierr
;
end
