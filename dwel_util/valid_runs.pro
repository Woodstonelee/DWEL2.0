function valid_runs,valid,start,num
  compile_opt idl2
  ;
  ;valid_runs defines the separated "runs" of adjacent points from a logocal array
  ;of valid points. It is used to fill invalid points in DWEL processing but is a general routine.
  ;     valid is a logical array with 1 as valid and 0 as not
  ;     start is output and has the starting positions of the runs in valid
  ;     num is out put. OIt is thae same length as start and has the number of points in the run
  ;
  ;     If the function returns "false" it means there are no valid cases at all
  ;
  start=[-1]
  num=[-1]
  pos_valid=where(valid,npos_valid)
  if (npos_valid le 0) then return,0b
  delta=pos_valid-shift(pos_valid,+1)
  ;print,'valid=['+strjoin(strtrim(valid,2),',',/single)+']'
  ;print,'pos_valid=['+strjoin(strtrim(string(pos_valid),2),',',/single)+']'
  ;print,'delta=['+strjoin(strtrim(string(delta),2),',',/single)+']'
  ;initialize the accumulators
  num=[1]
  j=0
  start=[pos_valid[0]]
  if (npos_valid gt 1) then begin
    for k=1,npos_valid-1 do begin
      if (delta[k] eq 1) then begin
        num[j]=num[j]+1
      endif else begin
        num=[num,1]
        j=j+1
        start=[start,pos_valid[k]]
      endelse
    endfor
  endif
  ;print,'start=['+strjoin(strtrim(string(start),2),',',/single)+']'
  ;print,'num=['+strjoin(strtrim(string(num),2),',',/single)+']'
  return,1b
end
