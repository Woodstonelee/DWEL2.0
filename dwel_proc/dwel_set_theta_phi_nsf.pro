function dwel_set_theta_phi_nsf, pstat,zen_tweak
  compile_opt idl2
  
  resolve_routine, 'TAG_EXIST', /compile_full_file, /either

  ;; see if the user has provided a scan encoder of zenith point
  if tag_exist((*pstat), 'ZenEnc') then begin
    zenith_count = (*pstat).ZenEnc + long(zen_tweak)
  endif else begin
    zenith_count = 392848L + long(zen_tweak) ; NSF DWEL
  endelse 
  
  ShotZen=fltarr((*pstat).NShots,(*pstat).Nscans)
  ShotAzim=fltarr((*pstat).Nshots,(*pstat).Nscans)
  
  ;compute the zenith [-180,180] and azimuth [0,360]
  ;for every shot from the encoder values
  ;ShotZen = (262144.0-(*pstat).ShotZen)*360.0/524288.0
  ;ShotAzim = (*pstat).ShotAzim*360.0/524288.0
  ShotZen = ((double(zenith_count) - double((*pstat).ShotZen)) / double(524288)) * 360.0d0
  ShotAzim = ((double(524288)-double((*pstat).ShotAzim))/double(524288)) * 360.0d0
  
  ;find the wrap around due to offset
  index=where(ShotZen lt -180.0,count)
  if (count gt 0) then ShotZen[index]=ShotZen[index]+360.0d0
  index=0b
  index = where(ShotZen gt 180.0, count)
  if (count gt 0) then ShotZen[index] = ShotZen[index] - 360.0d0
  index = 0b
  ;fix convention of zenith and azimuth
  index = where(ShotZen lt 0.0, count)
  if (count gt 0) then begin
    ShotZen[index] = - ShotZen[index]
    ShotAzim[index] = ShotAzim[index] + 180.0d0
  endif
  index=0b
  ;shift azimuth to [0,360]
  index = where(ShotAzim[*] gt 360.0, count)
  if count gt 0 then begin
    ShotAzim[index] = ShotAzim[index] - 360.0
  endif
  index=0b
  ;Put the results back into the structure
  (*pstat).ShotZen=float(ShotZen)
  (*pstat).shotAzim=float(ShotAzim)
  
  ShotZen=0b
  ShotAzim=0b
  
  return,1b
  
end

