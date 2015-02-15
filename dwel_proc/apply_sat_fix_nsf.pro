function apply_sat_fix_nsf, basefixed_satwf, pulse_model, i_val, scale, satfixedwf=satfixedwf
  ; set some values - i_val now read in
  p_centpeak=i_val[2]
  p_troughloc=i_val[3]
  p_scdpeakloc=i_val[4]
  
  ; initial maximum peak location
  tmp = max(basefixed_satwf, maxpeakloc)
  allmaxpeakloc = where(basefixed_satwf eq tmp)
  maxpeakloc = maxpeakloc[0]
  nu_pos=where(basefixed_satwf eq tmp,nu_npos)
  if (nu_npos le 1) then return,2
  if (nu_npos gt 20) then return,0
  nu_pos=0b
  ; find the three zero-cross points after the maximum peak
  wflen = size(basefixed_satwf, /n_elements)
  zero_xloc = $
    where(basefixed_satwf[maxpeakloc:wflen-2]*basefixed_satwf[maxpeakloc+1:wflen-1] $
    le 0, tmpcount) + maxpeakloc
    
  if (size(zero_xloc, /n_elements) lt 3) then begin
    return, 0
  endif
  ; find the minimum and maximum between the first zero-cross point and the third zero-cross point
  tmp = min(basefixed_satwf[zero_xloc[0]:zero_xloc[1]],tmploc)
  ;check if the minimum has saturated at original zero DN
  ;if so, use mean of positions
  nu_pos=where(basefixed_satwf[zero_xloc[0]:zero_xloc[1]] eq tmp,nu_npos)
  if (nu_npos gt 1) then begin
    troughloc = round(mean(nu_pos,/double)) + zero_xloc[0]
  endif else troughloc=tmploc+zero_xloc[0]
  ; now get the max value in the second interval (note use of separate intervals)
  tmp = max(basefixed_satwf[zero_xloc[1]:zero_xloc[2]], tmploc)
  scdpeakloc = fix(mean(tmploc[0])) + zero_xloc[1]
  ; satfix scale now more stable
  satfix_scale = (2.0*abs(basefixed_satwf[troughloc[0]])+abs(basefixed_satwf[scdpeakloc[0]]))/ $
    (2.0*abs(pulse_model[p_troughloc])+abs(pulse_model[p_scdpeakloc]))
  ;
  ;now find the centre position if there is a run of top values of the main saturated waveform
  ;again use mean of positions as default position of the main peak
  ;use the mean of saturated positions
  centpeak = mean(allmaxpeakloc)

  ;now only use the 1% and 99% points rather than overwriting full pulse length
  plen = n_elements(pulse_model[i_val[1]:i_val[5]])
  satfixedwf = basefixed_satwf
  ;
  ;stable estimate of start uses estimates of main peak and other two with weights
  nl=round((3.0*float(centpeak-p_centpeak)+2.0*float(troughloc-p_troughloc) + $
    float(scdpeakloc-p_scdpeakloc))/6.0) + i_val[1]
  ;nr=troughloc+plen-p_troughloc-1
  nr=nl+plen-1
  nsatfix=n_elements(satfixedwf)
  if ((nl lt 0) or (nr gt nsatfix-1)) then return,0
  satfixedwf[nl:nr] = pulse_model[i_val[1]:i_val[5]]*satfix_scale
  temp=reform(satfixedwf[nl:nr])
  ;prevent overflow values if possible
  nu_pos=where(abs(temp) gt 32767.0/scale,nu_npos)
  if (nu_npos gt 0) then begin
    temp[nu_pos]=32767.0/scale
    satfixedwf[nl:nr]=temp
    return,0
  endif
  temp=0b
  nu_pos=0b
  return, 1
end
