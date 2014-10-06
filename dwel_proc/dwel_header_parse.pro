function DWEL_header_parse, info
compile_opt idl2
;
ipoint=0
if (n_elements(info) le 0) then return,ipoint
kk=0
accum=''
nn=n_elements(info)
for i=0,nn-1 do begin
  info[i]=strtrim(info[i],2)
  np=strpos(info[i],'=')
  if(kk le 0) then begin
    if (np lt 0) then begin
      info[ipoint]=info[i]
      ipoint=ipoint+1
      goto,endloop
    endif else begin
      buf=strtrim(strmid(info[i],np+1),2)
      if (strpos(buf,'(') eq 0) then begin
        if (strpos(buf,')') ne strlen(buf)-1) then begin
          kk=kk+1
          accum=strtrim(info[i],2)
          goto,endloop
        endif
      endif else begin
        info[ipoint]=info[i]
        ipoint=ipoint+1
        goto,endloop
      endelse
    endelse
  endif else begin
    if (strpos(info[i],')') lt 0) then begin
      accum=accum+','+strtrim(info[i],2)
      kk=kk+1
      goto,endloop
    endif else begin
      info[ipoint]=accum+','+strtrim(info[i],2)
      kk=0
      accum=''
      ipoint=ipoint+1
      goto,endloop
    endelse
  endelse
endloop:
endfor
info=info[0:ipoint-1]
;
return,ipoint
;
end
