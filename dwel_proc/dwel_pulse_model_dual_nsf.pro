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
      
    pulse=[ $
      0.003419, 0.003391, 0.003364, 0.003337, 0.003310, 0.003282, -0.000085, $
      -0.000113, -0.000140, -0.003507, -0.003535, -0.003562, -0.003589, -0.000276, $
      -0.003644, -0.003671, -0.000358, -0.000385, -0.000412, 0.002901, 0.002874, $
      0.002846, 0.002819, 0.002792, 0.002765, 0.002737, 0.002710, 0.002683, $
      -0.000685, -0.000712, -0.000739, -0.000766, -0.000794, -0.000821, -0.000848, $
      -0.000875, -0.000903, -0.000930, -0.000957, -0.000984, 0.002329, 0.002302, $
      0.002274, 0.002247, 0.002220, -0.001148, -0.001175, -0.001202, -0.001230, $
      -0.004597, -0.001284, -0.001311, -0.004679, -0.001366, -0.001393, 0.001920, $
      0.001893, 0.001866, 0.001838, 0.001811, 0.001784, -0.001584, -0.001611, $
      -0.001638, -0.001665, -0.005033, -0.005060, -0.005088, -0.005115, 0.001539, $
      0.008192, 0.024867, 0.071604, 0.158426, 0.295352, 0.448980, 0.632672, $
      0.792981, 0.923227, 0.993347, 1.000000, 0.986611, 0.889714, 0.736032, $
      0.538924, 0.338477, 0.154731, 0.004389, -0.125912, -0.222809, -0.259580, $
      -0.266288, -0.249613, -0.199535, -0.129416, -0.059296, 0.000803, 0.050881, $
      0.084257, 0.097591, 0.104245, 0.104217, 0.087488, 0.064079, 0.040669, $
      0.013919, -0.002810, -0.019539, -0.029587, -0.036295, -0.036322, -0.029669, $
      -0.023015, -0.019702, -0.006368, -0.003055, 0.006939, 0.010252, 0.010225, $
      0.013538, 0.010170, 0.010143, 0.010116, 0.006748, 0.003381, 0.003353, $
      -0.000014, -0.000041, -0.003409, -0.003436, -0.006804, -0.006831, -0.006858, $
      -0.006886, -0.003573, -0.006940, -0.003627, -0.000314, -0.000341, 0.002972, $
      0.006285, 0.002917, 0.006231, 0.006203, 0.006176, 0.002808, -0.000559, $
      -0.000586, -0.003954, -0.003981, -0.007349, -0.010716, -0.007403, -0.007430, $
      -0.007458, -0.004145, -0.000832, 0.005822, 0.005795, 0.009108, 0.009081 $
      ]
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
      0.003208, 0.003187, -0.003536, 0.003144, 0.003123, 0.003101, 0.003080, $
      0.003058, -0.003664, 0.003016, -0.003707, -0.003729, -0.003750, 0.002930, $
      -0.003793, 0.002887, -0.003835, -0.003857, 0.002823, -0.003900, -0.003921, $
      0.002759, 0.002738, -0.003985, -0.004007, -0.004028, -0.004049, -0.004071, $
      -0.004092, -0.004113, 0.002566, 0.009246, 0.002524, 0.002502, 0.002481, $
      0.002460, 0.002438, 0.002417, -0.004306, 0.002374, -0.004349, -0.011072, $
      -0.004392, -0.011114, -0.004434, -0.004456, 0.002224, 0.002203, 0.002182, $
      0.002160, 0.008840, 0.008819, 0.008797, 0.002075, 0.002053, 0.002032, $
      -0.004691, -0.004712, 0.001968, -0.004755, -0.004777, -0.004798, 0.001882, $
      0.001861, 0.001839, 0.001818, 0.001797, 0.001775, 0.001754, 0.001732, $
      0.001711, 0.008391, 0.021772, 0.055258, 0.142354, 0.289763, 0.470678, $
      0.658295, 0.839210, 0.959813, 1.000000, 0.979875, 0.912840, 0.758687, $
      0.530819, 0.262743, 0.021473, -0.152784, -0.313638, -0.387374, -0.394097, $
      -0.380716, -0.313724, -0.219926, -0.112725, -0.012226, 0.068169, 0.121758, $
      0.155243, 0.155222, 0.148499, 0.128374, 0.094846, 0.061317, 0.021088, $
      -0.012440, -0.032566, -0.045990, -0.052713, -0.052734, -0.052755, -0.039374, $
      -0.025993, -0.005910, 0.007471, 0.014151, 0.027533, 0.034213, 0.034191, $
      0.027468, 0.020746, 0.020724, 0.007300, 0.000577, -0.006145, -0.012868, $
      -0.019591, -0.012911, -0.012932, -0.006252, -0.012975, 0.000406, 0.007086, $
      0.007065, 0.007044, 0.007022, 0.007001, 0.006979, 0.006958, 0.006937, $
      0.006915, 0.000192, 0.000171, -0.006552, -0.006573, -0.006594, -0.006616, $
      -0.006637, 0.000043, -0.006680, 0.000000, -0.000021, -0.000043, 0.006637, $
      0.006616, 0.006594, 0.006573, -0.000150, -0.000171, 0.006509, 0.006487 $
      ]
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

