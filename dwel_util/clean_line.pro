pro clean_line, inline,val_point,outline,filt,err
  compile_opt idl2
  ;clean line cleans a line by smoothing with filt
  ;the null places are first filled by interpolation or extrapolation
  ;the length of filt determines the length of an initial median smooth
  ;finally it is smoothed by filt
  ;
  ;inline is the line of data
  ;val_point is a logical indicator for the valid points
  ;outline is the result (same length as inline)
  ;filt is the filter (eg gaussian)
  ;err flags an error
  ;

  resolve_routine, 'VALID_RUNS', /compile_full_file, /either

  ;start with setup and consistency
  err=0
  npoint=n_elements(inline)
  if (npoint le 1) then begin
    err=1
    goto,out
  endif
  flen=n_elements(filt)
  if (flen le 2) then begin
    err=2
    goto,out
  endif
  frad=(flen-1)/2
  if (n_elements(val_point) ne npoint) then begin
    err=3
    goto,out
  endif
  outline=fltarr(npoint)
  
  ;fill the holes in the line
  ;first get the runs of places to fill
  istat=valid_runs(~val_point,start,num)
  if (istat) then begin
    ;now for the interpolation
    num_runs=n_elements(start)
    outline=float(inline)
    num_vec=n_elements(outline)
    for k=0,num_runs-1 do begin
      if (start[k] eq 0) then begin
        outline[start[k]:start[k]+num[k]-1]=outline[start[k]+num[k]]
      endif else if ((start[k]+num[k]) ge num_vec) then begin
        outline[start[k]:start[k]+num[k]-1]=outline[start[k]-1]
      endif else begin
        int_pos=indgen(num[k])
        places=replicate(1.0,num[k])-(float(int_pos)+replicate(1.0,num[k]))/float(num[k]+1)
        outline[start[k]:start[k]+num[k]-1]=outline[start[k]-1]*places+ $
          outline[start[k]+num[k]]*(replicate(1.0,num[k])-places)
      endelse
    endfor
  endif else begin
    outline=float(inline)
  endelse

  ;; temp=[replicate(outline[0],2*frad),float(outline),replicate(outline[npoint-1],2*frad)]
  ;; temp=median(temp,flen,/double)
  ;; temp=convol(temp,filt)
  ;; outline=temp[2*frad:npoint+2*frad-1]

  ;
  ;; flen is usually too small for median filter to remove noise. There could be
  ;; five adjacent bins of abnormal values in DWEL (NSF) data. A flen of 20 is
  ;; usually good. While flen is currently chosen as 5, we set the width of
  ;; median filter as flen*4
  medfwidth = flen*4
  temp=[replicate(outline[0],medfwidth),float(outline),replicate(outline[npoint-1],medfwidth)]
  temp=median(temp,medfwidth,/double)
  temp = temp[medfwidth:npoint+medfwidth-1]
  ;; to avoid really weird value in the first few lines or last few lines
  temp[0:fix(medfwidth)/2-1] = temp[fix(medfwidth/2)]
  temp[npoint-fix(medfwidth/2):npoint-1] = temp[npoint-fix(medfwidth/2)-1]
  temp=[replicate(temp[0],2*frad),float(temp),replicate(temp[npoint-1],2*frad)]
  temp=convol(temp,filt)
  outline=temp[2*frad:npoint+2*frad-1]
  temp=0b
  out:
  return
end
