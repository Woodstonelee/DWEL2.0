function dwel_set_theta_phi_oz, pstat,zen_tweak
compile_opt idl2

;zenith_count=185896L ;old guess by dlbj
;zenith_count=185600L ;best early (June)
;zenith_count=184800L ;maybe best overall to end of July
;zenith_count=182196L ; new estimate 01082014 by MS after changes
zenith_count=182196L+long(zen_tweak)

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

