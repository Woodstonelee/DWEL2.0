pro DWEL_pulse_model_dual_oz, wavelength, i_val, t_val, r_val, p_range, p_time, pulse, t_fwhm, r_fwhm
compile_opt idl2

;THis routine sets up the model for the base recorded DWEL pulse. 
;Pulses have been set as an empirical model as an array at the step corresponding to 0.5 nsec sampling.
;This model may be interpolated to fit any other step size and is controlled by the position and magnitude
;of the peak of the pulse - the T_zero point
;
;Some information:
;
; Peak occurs at zero range or time with a value of 1.0
; FWHM is the sum of pulse values
;
;INPUTS:
;
;     wavelength  wavelength of the laser of which the pulse model is needed
;     i_val     positions of first,Left 0.0001,peak,trough,secondary peak,Right 0.0001, last points
;     t_val     time vale of first,Left 0.0001,peak,trough,secondary peak,Right 0.0001, last points
;     r_val     range vale of first,Left 0.0001,peak,trough,secondary peak,Right 0.0001, last points
;     p_range   the range ordinates of the model
;     p_time    the time ordinates of the model
;     pulse     values of the normalised pulse (peak value = 1.0)
;

c=0.299792458
c2=c/2.0

;; Oz DWEL July 2014
if (wavelength eq 1064) then begin
; Mean pulse model
p_time=[$
-40.00,-39.50,-39.00,-38.50,-38.00,-37.50,-37.00,-36.50,-36.00,-35.50,$
-35.00,-34.50,-34.00,-33.50,-33.00,-32.50,-32.00,-31.50,-31.00,-30.50,$
-30.00,-29.50,-29.00,-28.50,-28.00,-27.50,-27.00,-26.50,-26.00,-25.50,$
-25.00,-24.50,-24.00,-23.50,-23.00,-22.50,-22.00,-21.50,-21.00,-20.50,$
-20.00,-19.50,-19.00,-18.50,-18.00,-17.50,-17.00,-16.50,-16.00,-15.50,$
-15.00,-14.50,-14.00,-13.50,-13.00,-12.50,-12.00,-11.50,-11.00,-10.50,$
-10.00,-9.50,-9.00,-8.50,-8.00,-7.50,-7.00,-6.50,-6.00,-5.50,$
-5.00,-4.50,-4.00,-3.50,-3.00,-2.50,-2.00,-1.50,-1.00,-0.50,$
0.00,0.50,1.00,1.50,2.00,2.50,3.00,3.50,4.00,4.50,$
5.00,5.50,6.00,6.50,7.00,7.50,8.00,8.50,9.00,9.50,$
10.00,10.50,11.00,11.50,12.00,12.50,13.00,13.50,14.00,14.50,$
15.00,15.50,16.00,16.50,17.00,17.50,18.00,18.50,19.00,19.50,$
20.00,20.50,21.00,21.50,22.00,22.50,23.00,23.50,24.00,24.50,$
25.00,25.50,26.00,26.50,27.00,27.50,28.00,28.50,29.00,29.50,$
30.00,30.50,31.00,31.50,32.00,32.50,33.00,33.50,34.00,34.50,$
35.00,35.50,36.00,36.50,37.00,37.50,38.00,38.50,39.00,39.50,$
40.00]

pulse=[$
0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,-0.000676,-0.000676,-0.000676,$
-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,$
-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,$
-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,$
-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,$
-0.000676,-0.000676,-0.000676,-0.000676,-0.000676,-0.000661,-0.000663,-0.000677,-0.000755,-0.000629,$
-0.000584,-0.000455,-0.000244,0.000119,0.000294,0.000566,0.000879,0.003723,0.005502,0.011154,$
0.021713,0.044893,0.085870,0.154438,0.254019,0.387290,0.558426,0.722474,0.866768,0.967142,$
1.000000,0.956108,0.842129,0.675911,0.474706,0.265679,0.070073,-0.092647,-0.213085,-0.284504,$
-0.309195,-0.291783,-0.243520,-0.174886,-0.099198,-0.024835,0.038396,0.087039,0.117281,0.131816,$
0.130534,0.117319,0.094553,0.068864,0.039836,0.012677,-0.010774,-0.028513,-0.040564,-0.046049,$
-0.046414,-0.040637,-0.032580,-0.020497,-0.010042,0.000943,0.009570,0.015681,0.019465,0.020269,$
0.017930,0.014277,0.008879,0.005369,0.000879,-0.001622,-0.002992,-0.002905,-0.002392,-0.001328,$
-0.001085,0.000483,0.000851,0.002152,0.002095,0.002362,0.002378,0.002119,0.001670,0.001726,$
0.001450,0.001365,0.001048,0.001266,0.000995,0.001546,0.001781,0.001907,0.001830,0.001641,$
0.001456,0.001853,0.001393,0.001068,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,$
0.000000]
endif 
if (wavelength eq 1548) then begin
; Mean pulse model
p_time=[$
-40.00,-39.50,-39.00,-38.50,-38.00,-37.50,-37.00,-36.50,-36.00,-35.50,$
-35.00,-34.50,-34.00,-33.50,-33.00,-32.50,-32.00,-31.50,-31.00,-30.50,$
-30.00,-29.50,-29.00,-28.50,-28.00,-27.50,-27.00,-26.50,-26.00,-25.50,$
-25.00,-24.50,-24.00,-23.50,-23.00,-22.50,-22.00,-21.50,-21.00,-20.50,$
-20.00,-19.50,-19.00,-18.50,-18.00,-17.50,-17.00,-16.50,-16.00,-15.50,$
-15.00,-14.50,-14.00,-13.50,-13.00,-12.50,-12.00,-11.50,-11.00,-10.50,$
-10.00,-9.50,-9.00,-8.50,-8.00,-7.50,-7.00,-6.50,-6.00,-5.50,$
-5.00,-4.50,-4.00,-3.50,-3.00,-2.50,-2.00,-1.50,-1.00,-0.50,$
0.00,0.50,1.00,1.50,2.00,2.50,3.00,3.50,4.00,4.50,$
5.00,5.50,6.00,6.50,7.00,7.50,8.00,8.50,9.00,9.50,$
10.00,10.50,11.00,11.50,12.00,12.50,13.00,13.50,14.00,14.50,$
15.00,15.50,16.00,16.50,17.00,17.50,18.00,18.50,19.00,19.50,$
20.00,20.50,21.00,21.50,22.00,22.50,23.00,23.50,24.00,24.50,$
25.00,25.50,26.00,26.50,27.00,27.50,28.00,28.50,29.00,29.50,$
30.00,30.50,31.00,31.50,32.00,32.50,33.00,33.50,34.00,34.50,$
35.00,35.50,36.00,36.50,37.00,37.50,38.00,38.50,39.00,39.50,$
40.00]

pulse=[$
0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,-0.000686,-0.000686,-0.000686,$
-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,$
-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,$
-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,$
-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,-0.000686,$
-0.000686,-0.000686,-0.000686,-0.000686,-0.000685,-0.000698,-0.000635,-0.000731,-0.000596,-0.000805,$
-0.000733,-0.000927,-0.000758,-0.000925,-0.000707,-0.001090,-0.000531,-0.000197,0.001228,0.003674,$
0.011955,0.029858,0.067044,0.130567,0.230782,0.365109,0.528668,0.699771,0.856538,0.963581,$
1.000000,0.954043,0.834670,0.653477,0.439397,0.217075,0.015368,-0.149627,-0.263098,-0.322472,$
-0.327448,-0.290522,-0.222264,-0.140889,-0.057206,0.015152,0.073071,0.108630,0.127622,0.126521,$
0.114002,0.088731,0.062432,0.033568,0.009267,-0.012873,-0.026238,-0.036036,-0.038030,-0.037158,$
-0.029540,-0.022066,-0.011747,-0.004125,0.005592,0.011358,0.016526,0.018291,0.019501,0.017321,$
0.015636,0.010075,0.006862,0.002375,-0.000191,-0.003526,-0.003636,-0.004324,-0.002176,-0.002730,$
0.000228,0.000497,0.002135,0.002439,0.003993,0.004272,0.004562,0.003074,0.003585,0.002206,$
0.002659,0.000869,0.001604,0.000629,0.001822,0.000050,0.001825,0.001211,0.002769,0.001810,$
0.002806,0.002038,0.002163,0.000580,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,$
0.000000]
endif

p_range = p_time*c2
fwhm=total(pulse,/double)
t_fwhm=fwhm*(p_time[1]-p_time[0])    ; nanosec
r_fwhm=fwhm*(p_range[1]-p_range[0])  ; metres

;find where the left and right 0.0001*peak is
tmp_I = where(pulse ge 0.0001)
tmp = max(pulse, max_I)

;=====================================
; find the trough and secondary peak
; maximum peak location
tmp = max(pulse, maxpeakloc)
maxpeakloc = maxpeakloc[0]
; find the three zero-cross points after the maximum peak
wflen = size(pulse, /n_elements)
zero_xloc = where(pulse[maxpeakloc:wflen-2]*pulse[maxpeakloc+1:wflen-1] le 0, tmpcount) + maxpeakloc
; find the minimum and maximum between the first zero-cross point and the third zero-cross point
tmp = min(pulse[zero_xloc[0]:zero_xloc[2]], tmploc)
troughloc = fix(mean(tmploc)) + zero_xloc[0]
tmp = max(pulse[zero_xloc[0]:zero_xloc[2]], tmploc)
scdpeakloc = fix(mean(tmploc)) + zero_xloc[0]
;=====================================

i_val=[0, tmp_I[0],  max_I, troughloc, scdpeakloc, tmp_I[n_elements(tmp_I)-1], n_elements(pulse)-1]

r_val=[p_range[i_val[0]],p_range[i_val[1]],p_range[i_val[2]], p_range[i_val[3]], p_range[i_val[4]], p_range[i_val[5]],p_range[i_val[6]]]
t_val=r_val/c2

print,'p_val check=',[pulse[i_val[0]],pulse[i_val[1]],pulse[i_val[2]],pulse[i_val[3]],pulse[i_val[4]],pulse[i_val[5]],pulse[i_val[6]]]
print,'expected number=',i_val[6]+1,' actual number=',n_elements(pulse)

return
end

