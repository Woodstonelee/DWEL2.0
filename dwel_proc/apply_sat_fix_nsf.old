function apply_sat_fix_nsf, basefixed_satwf, pulse_model, i_val, scale_mean, satfixedwf=satfixedwf
  compile_opt idl2
  ; set some values - i_val now read in
  p_centpeak=i_val[2]
  p_troughloc=i_val[3]
  p_scdpeakloc=i_val[4]
  ;
  p_index=indgen(n_elements(pulse_model))
  ; initial maximum peak location
  tmp = max(basefixed_satwf, maxpeakloc)
  maxpeakloc = maxpeakloc[0]
  hi_pos=where(basefixed_satwf ge tmp,hi_npos)
  if (hi_npos le 1) then return,2
  ;  if ((hi_npos gt 30) and ((hi_pos[hi_npos-1]-hi_pos[0]+1) gt 161)) then return,0
  if ((hi_pos[hi_npos-1]-hi_pos[0]) gt hi_npos+2) then return,0
  maxpeakloc=round(hi_pos[0])
  maxpeakmid=round(mean(hi_pos,/double))
  ; find the three zero-cross points after the maximum peak
  wflen = size(basefixed_satwf, /n_elements)
  zero_xloc = where(basefixed_satwf[maxpeakloc:wflen-2]*basefixed_satwf[maxpeakloc+1:wflen-1] le 0, tmpcount) + maxpeakloc
  if (size(zero_xloc, /n_elements) lt 3) then begin
    return, 0
  endif
  ; find the minimum and maximum between the first zero-cross point and the third zero-cross point
  tmp = min(basefixed_satwf[zero_xloc[0]:zero_xloc[1]],tmploc)
  ;check if the minimum has saturated at original zero DN
  ;if so, use mean of positions
  ;lo_pos is the list of saturate negative lobe points if lo_npos gt 1
  lo_pos=where(basefixed_satwf[zero_xloc[0]:zero_xloc[1]] eq tmp,lo_npos)
  if (lo_npos gt 1) then begin
    troughloc = round(mean(lo_pos,/double)) + zero_xloc[0]
  endif else troughloc=tmploc+zero_xloc[0]
  lo_LH=lo_pos[0]+zero_xloc[0]-1
  lo_RH=lo_pos[lo_npos-1]+zero_xloc[0]+1
  p_troughloc_LH=p_troughloc+(lo_LH-troughloc)
  p_troughloc_RH=p_troughloc+(lo_RH-troughloc)
  
  ; now get the max value in the second interval (note use of separate intervals)
  tmp = max(basefixed_satwf[zero_xloc[1]:zero_xloc[2]], tmploc)
  scdpeakloc = fix(mean(tmploc[0])) + zero_xloc[1]
  ;
  ;now find the main peak centre position if there is a run of top values of the main saturated waveform
  ;again use mean of positions as default position of the main peak
  tmp = max(basefixed_satwf[maxpeakloc:zero_xloc[0]],cloc)
  hi_pos=where(basefixed_satwf[maxpeakloc:zero_xloc[0]] ge tmp,hi_npos)
  if (hi_npos gt 1) then begin
    centpeak=mean(float(hi_pos),/double)+float(maxpeakloc)
  endif else centpeak=float(maxpeakloc+cloc)
  centloc=round(centpeak)
  hi_LH=hi_pos[0]+maxpeakloc-1
  hi_RH=hi_pos[hi_npos-1]+maxpeakloc+1
  ;
  ;stable estimate of start uses estimates of main peak and other two with weights
  fpos=(3.0*float(centpeak-p_centpeak)+2.0*float(troughloc-p_troughloc)+ $
    float(scdpeakloc-p_scdpeakloc))/6.0
  ;
  ;now resample the pulse so that the peak is at the right place
  ;======================================================
  ;new code
  npulse=n_elements(pulse_model)
  wp_out=fltarr(npulse)
  x_in=fpos+float(p_index)
  x_out=float(p_index+round(fpos))
  wp_out=interpol(pulse_model,x_in,x_out)
  x_out=round(x_out)
  x_in=0b
  ;end of nw code
  ;======================================================
  
  p_centpeak_test=round(float(centpeak)-fpos)
  tmp=max(wp_out,p_centpeak_nu)
  p_hi_LH=p_centpeak_nu+round(hi_LH-centloc)
  p_hi_RH=p_centpeak_nu+round(hi_RH-centloc)
  
  ;satfix scale is more stable after changes!
  satfix_scale = $
    ((abs(basefixed_satwf[hi_LH])+abs(basefixed_satwf[hi_RH]))+abs(basefixed_satwf[lo_LH])+abs(basefixed_satwf[lo_RH])+abs(basefixed_satwf[scdpeakloc]))/ $
    ((abs(wp_out[p_hi_LH])+abs(wp_out[p_hi_RH]))+abs(wp_out[p_troughloc_LH])+abs(wp_out[p_troughloc_RH])+abs(wp_out[p_scdpeakloc]))
    
  ;now only use the 1% and 99% points rather than overwriting full pulse length
  plen = n_elements(pulse_model[i_val[1]:i_val[5]])
  satfixedwf = basefixed_satwf
  
  nl=x_out[0]+i_val[1]
  nr=nl+plen-1
  nsatfix=n_elements(satfixedwf)
  if ((nl lt 0) or (nr gt nsatfix-1)) then begin
    return,0
  endif
  satfixedwf[nl:nr] = wp_out[i_val[1]:i_val[5]]*satfix_scale
  temp=reform(satfixedwf[nl:nr])
  ;prevent overflow values if possible
  nu_pos=where(abs(temp) gt 32767.0/scale_mean,nu_npos)
  if (nu_npos gt 0) then begin
    temp[nu_pos]=32767.0/scale_mean
    satfixedwf[nl:nr]=temp
    return,0
  endif
  temp=0b
  nu_pos=0b
  wp_out=0b
  return, 1
end
