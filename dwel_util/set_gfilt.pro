function set_gfilt,filter_rad
compile_opt idl2
;
filter_len=2L*long(filter_rad)+1L
filt=fltarr(filter_len)
sig=0.329505*float(filter_rad)
for j=0L,filter_len-1L do begin
  filt[j]=exp(-0.5*(float(j-filter_rad)/sig)^2)
endfor
filt=filt/total(filt)
return,filt
;
end
