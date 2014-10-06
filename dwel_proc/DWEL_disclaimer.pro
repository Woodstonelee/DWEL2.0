function DWEL_disclaimer, rev
compile_opt idl2

;The disclaimer hopefully makes claims to the good bits
;but absolves of any blame due to bad bits

text=[ $
'Copyright(C) CSIRO(Aust) & Boston Uni, 2012-2014, V1 Rev'$
+strtrim(string(rev,format='(i9)'),2),$
'Provided by CSIRO Marine & Atmospheric Sciences and Boston',$
'University extending original previous CSIRO copyrights',$
'Copyright(C) CSIRO CMAR Canopy Lidar Initiative 2004-2012,',$
'Copyright(C) CSIRO Earth Observation Centre 2002-2004.',$
'This Module is provided for the use of sophisticated',$
'users who have participated in a training course',$
'or are collaborating with CSIRO and BU in joint research projects',$
'and who have the means to make independent assessments',$
'of the suitability of the algorithms provided.',$
'Due to limited internal review due to the developmental',$
'nature of the work, recommendations, equations, algorithms,',$
'calculations and/or implementations in this software are',$
'provided only for the purposes of research, teaching',$
'and information, and should not be taken as final',$
'publication or commercial quality and validity.',$
'CSIRO and BU do not offer any support or guarantee for this code',$
'outside of its use for Joint research involving them or in training courses.' $
]

return,text

end
