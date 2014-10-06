;======================================================================
;  The following set of procedures is part of the CSIRO EOC Hyperion
;  Workshop module set. The routines all have a similar structure.
;  They can be used stand-alone or as part of a Project structure with
;  ENVI main menu button.
;
;  The general structure for a procedure called "name" is:
;
;  [Procedures called from name-doit]
;    These may be called anything
;    they arerranged so that every procedure occurs before its first
;    reference. This is good IDL practice.
;  pro name_doit
;    the main area of action for the algorithm - the DOIT!
;  [Procedures supporting the general top level procedure]
;      These procedures are not part of the algorithm - just the framework.
;      The workshop routines have a common top level look and feel.
;    pro name_banner_text
;      Top level widget Help and how to use the module
;    pro name_disclaimer
;      General disclaimer to all responsibility and blame accessed
;      from main widget
;    pro name_action_go
;      Main widget event handler
;    pro name_action_exit
;      Another event handler
;    pro name_resize
;      Yet another event handler
;    pro name_cleanup
;      Will they never end?
;    pro name
;      This is the common top level procedure that simply sets up an
;      IDL widget structure that allows the user to look at the main
;      help and if desired hit GO and start the module. It handles
;      exit nicely.
;======================================================================
;

function set_ptcl_base_params, pstat
compile_opt idl2

set_params:

;first set up the widget base and user help
wb_fix_ifov=widget_auto_base( $
  title='Set parameters for Point Cloud Sieve', $
  group=(*pstat).top)

Info_Text=[$
'The following information allows you to sieve the points cloud data',$
'on x, y, z and intensity. It asks for simple ranges in these parameters',$
'and the program then sieves the points satisfying these ranges from the',$
'input file.',$
' ' $
]

;put in the help and some next level bases to use
;for input info
w_fix_txt=widget_slabel(wb_fix_ifov,prompt=Info_Text,$
  /frame,xsize=43,ysize=5)

wb_lower=widget_base(wb_fix_ifov,/row,/frame)

wb_fix_par=widget_base(wb_lower,/row,/frame)
wb_fix_par_1=widget_base(wb_lower,/row,/frame)

;set the x direction
w_x_z=widget_param(wb_fix_par,/auto_manage,default=(*pstat).min_X,dt=type,$
         Prompt='Min_X: ', $
         uvalue='min_x',floor=-1000.0,ceil=1000.0,field=2)
w_y_z=widget_param(wb_fix_par,/auto_manage,default=(*pstat).max_X,dt=type,$
         Prompt='Max_X: ', $
         uvalue='max_x',floor=-1000.0,ceil=1000.0,field=2)

;set the y direction
w_x_R=widget_param(wb_fix_par_1,/auto_manage,default=(*pstat).min_Y,dt=type,$
         Prompt='Min_Y: ', $
         uvalue='min_y',floor=-1000.0,ceil=1000.0,field=2)
w_y_R=widget_param(wb_fix_par_1,/auto_manage,default=(*pstat).max_Y,dt=type,$
         Prompt='Max_Y: ', $
         uvalue='max_y',floor=-1000.0,ceil=1000.0,field=2)

;-----------------------------------

wb_bot=widget_base(wb_fix_ifov,/row,/frame)

wb_fix_par_3=widget_base(wb_bot,/row,/frame)
wb_fix_par_4=widget_base(wb_bot,/row,/frame)

;set the z direction
w_x_z=widget_param(wb_fix_par_3,/auto_manage,default=(*pstat).min_Z,dt=type,$
         Prompt='Min_Z: ', $
         uvalue='min_z',floor=-1000.0,ceil=1000.0,field=2)
w_y_z=widget_param(wb_fix_par_3,/auto_manage,default=(*pstat).max_Z,dt=type,$
         Prompt='Max_Z: ', $
         uvalue='max_z',floor=-1000.0,ceil=1000.0,field=2)

;set the I direction
w_x_R=widget_param(wb_fix_par_4,/auto_manage,default=(*pstat).min_I,dt=type,$
         Prompt='Min_I: ', $
         uvalue='min_i',floor=-1000.0,ceil=1000.0,field=2)
w_y_R=widget_param(wb_fix_par_4,/auto_manage,default=(*pstat).max_I,dt=type,$
         Prompt='Max_I: ', $
         uvalue='max_i',floor=0.0,ceil=30000.0,field=1)

;realise the widget using ENVI auto-manage
result=auto_wid_mng(wb_fix_ifov)

;result is a structure - see if setup cancelled
if (result.accept eq 0) then return,0b

if($
(result.min_x gt result.max_x) or $
(result.min_x gt result.max_x) or $
(result.min_x gt result.max_x) or $
(result.min_x gt result.max_x)) then begin
  info_text=['Inadmissible input (min must be <= max)',$
  'try again (Y/N) ?']
  result=dialog_message(info_text,/question,/default_no, $
  title='Bad input to the sieve !',$
  dialog_parent=(*pstat).top)
  if(result eq 'No') then begin
    return,0b
  endif else goto, set_params
endif

(*pstat).min_x=result.min_x
(*pstat).max_x=result.max_x
(*pstat).min_y=result.min_y
(*pstat).max_y=result.max_y
(*pstat).min_z=result.min_z
(*pstat).max_z=result.max_z
(*pstat).min_i=result.min_i
(*pstat).max_i=result.max_i

return,1b

end

pro pt_util_doit, event
compile_opt idl2

;The logic for pt_util is to:
;     Select a point cloud ascii file
;     Set sieve bounds
;     Sieve the and write out new LOG file (only)
;
;All user interaction is by widget

tfile=99

;First setup protective aunty-Catch to watch out for errors
error_status=0
catch, error_status
if (error_status ne 0) then begin
   catch,/cancel
   help,/last,output=out
   info_text=[strtrim(string('Catch trapped an Error !'),2),$
              strtrim('Error Name : '+strtrim(!err_string,2),2),$
              strtrim(string('Error Number: ',error_status,$
              format='(a,i5)'),2),$
              strtrim('Last Message: '+strtrim(out,2),2)]
   result=dialog_message(info_text,/error,title='Error in pt_util')
   goto, cleanup
  ;now see if the output file has a handle
  if(outfile ne '') then begin
    fids=envi_get_file_ids()
    if(fids[0] ne -1) then begin
      for i=0,n_elements(fids)-1 do begin
        envi_file_query,fids[i],fname=tname
        if (strtrim(outfile,2) eq $
          strtrim(tname,2)) then begin
          envi_file_mng,id=fids[i],/remove
        endif
      endfor
    endif
  endif
  free_lun,tfile,/force
endif

;get the grandfather widget in case...
widget_control,event.top,get_uvalue=pstate

;clean up any fids which are no longer where they were!
;ENVI issue that is annoying and leads to confusion
clean_envi_file_fids

;define the input and output images as strings
cal_name=''
out_name=''
o_name=''
tfile=99

cfg = envi_get_configuration_values()
if (cfg.DEFAULT_DATA_DIRECTORY eq '') then begin
  f_base='DWEL_Point_Cloud_File.log'
  f_path=strtrim(file_dirname(f_base),2)
endif else begin
  f_base='DWEL_Point_Cloud_File.log'
  f_path=strtrim(cfg.DEFAULT_DATA_DIRECTORY,2)
endelse

input:

;---------------------------------------------------------------------
;get a name for the ascii text file (must be a .txt file)
gain_file:

n_base=strlen(f_base)
n_dot=strpos(f_base,'.',/reverse_search)

if (cal_name eq '') then begin
  if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
    cal_name=f_base
  endif else begin
    cal_name=strmid(f_base,0,n_dot)+'.log'
  endelse
endif

cal_name=dialog_pickfile(title='Select DWEL ASCII Point Cloud (.log) file',$
         filter='*.log',file=cal_name, $
         dialog_parent=event.top, $
         path=f_path)
if(cal_name eq '') then begin
  result=dialog_message('Try again ? (No/Yes) ',/question,/default_no, $
  title='No point cloud file selected or selection cancelled !',$
  dialog_parent=event.top)
  if(result eq 'No') then begin
    goto, cleanup
  endif else goto, gain_file
endif

if(not file_test(cal_name)) then begin
  Info_Text=[$
  'There is no file called:',$
  strtrim(cal_name,2),$
  ' ',$
  'Try again ? (No/Yes) '$
  ]
  result=dialog_message(Info_Text,/question,/default_no, $
  title='Nominated point cloud file does not exist !',$
  dialog_parent=event.top)
  if(result eq 'No') then begin
    goto, cleanup
  endif else goto, gain_file
endif

f_base=cal_name

widget_control,/hourglass

;read the data
err=0
headers=['']
envi_read_cols, cal_name, temp, error=err,/read_skip,skip=headers

widget_control

if (err ne 0) then begin
  Info_Text=[$
  'There was a problem reading',$
  strtrim(cal_name,2),$
  ' ',$
  'Try another selection ? (No/Yes) '$
  ]
  result=dialog_message(Info_Text,/question,/default_no, $
  title='Error reading Base file',$
  dialog_parent=event.top)
  if(result eq 'No') then begin
    goto, cleanup
  endif else goto, gain_file
endif

nt=size(temp)

err=0
if((nt[0] gt 2) or (nt[0] le 1)) then begin
  Info_Text=[$
  'The Point Cloud file',$
  strtrim(cal_name,2),$
  'is NOT a 2D matrix',$
  ' ',$
  'Try a new file ? (No/Yes) '$
  ]
  result=dialog_message(Info_Text,/question,/default_no, $
  title='Incorrect Point Cloud file structure',$
  dialog_parent=event.top)
  if(result eq 'No') then begin
    goto, cleanup
  endif else goto, gain_file
endif

if(nt[1] lt 3) then begin
  Info_Text=[$
  'The Point Cloud file',$
  strtrim(cal_name,2),$
  'has too few columns - must at least have have (x,y,z)',$
  ' ',$
  'Try a new file ? (No/Yes) '$
  ]
  result=dialog_message(Info_Text,/question,/default_no, $
  title='Incorrect Point Cloud file structure',$
  dialog_parent=event.top)
  if(result eq 'No') then begin
    goto, cleanup
  endif else goto, gain_file
endif

;print,'matrix size=',nt[0],nt[1]
;print,'number of rows=',nt[2]
;print,'headers'
;for j=0,n_elements(headers)-1 do begin
;print,headers[j]
;endfor
;print,'first three rows'
;print,temp[*,0]
;print,temp[*,1]
;print,temp[*,2]

Num_Base=nt[2]
x=reform(temp[0,*])
y=reform(temp[1,*])
z=reform(temp[2,*])
intensity=reform(temp[3,*])

;print,'number of points in the input cloud=',n_elements(x)

n_points=n_elements(x)

Base_Info_Text=[ $
  'Point Cloud Data read from '+strtrim(cal_name,2),$
  ' ',$
  'Header(s) Read:',$
  headers, $
  ' ',$
  'Number of points in the Point Cloud='+strtrim(string(Num_Base,format='(i10)'),2),$
  'X Range (m)='+'['+strtrim(string(min(x),format='(f10.2)'),2)+','+strtrim(string(max(x),format='(f10.2)'),2)+']',$
  'Y Range (m)='+'['+strtrim(string(min(y),format='(f10.2)'),2)+','+strtrim(string(max(y),format='(f10.2)'),2)+']',$
  'Z Range (m)='+'['+strtrim(string(min(z),format='(f10.2)'),2)+','+strtrim(string(max(z),format='(f10.2)'),2)+']',$
  'I Range    ='+'['+strtrim(string(min(intensity),format='(f10.2)'),2)+','+strtrim(string(max(intensity),format='(f10.2)'),2)+']',$
  '',$
  'OK to go ahead (Y/N)?' $
  ]

min_x=min(x)
max_x=max(x)
min_y=min(y)
max_y=max(y)
min_z=min(z)
max_z=max(z)
min_i=min(intensity)
max_i=max(intensity)

result=dialog_message(Base_Info_Text,/question, $
title='Selected File OK to Sieve ?', dialog_parent=event.top)
if(result eq 'No') then begin
  goto,gain_file
endif

;-----------------------------------------------------------------------

;get parameters

;-----------------------------------------------------------------------

set_base:

;set up a structure and push it onto the heap
sav={ $
     min_x:min_x,$
     max_x:max_x,$
     min_y:min_y,$
     max_y:max_y,$
     min_z:min_z,$
     max_z:max_z,$
     min_i:min_i,$
     max_i:max_i,$
     top:event.top $
     }

;now locate the data on the heap with a pointer
p_stat=ptr_new(sav,/no_copy)

status=1b

;call to widget menu function
status = set_ptcl_base_params(p_stat)

if (not status) then begin
  ptr_free, p_stat
  result=dialog_message('No Parameters set in pt_util - Restart ? ',$
  /question, title='Parameter setup cancelled', $
  dialog_parent=event.top)
  if(result eq 'No') then begin
    goto, cleanup
  endif else goto, input
endif

;set the (valid) sizes back into the default
min_x=(*p_stat).min_x
max_x=(*p_stat).max_x
min_y=(*p_stat).min_y
max_y=(*p_stat).max_y
min_z=(*p_stat).min_z
max_z=(*p_stat).max_z
min_i=(*p_stat).min_i
max_i=(*p_stat).max_i
ptr_free, p_stat

output:

;get path and name as separate strings
last=strpos(cal_name,path_sep(),/reverse_search)
f_path=strmid(cal_name,0,last+1)
f_base=strtrim(strmid(cal_name,last+1,strlen(cal_name)-last-1),2)

n_base=strlen(f_base)
n_dot=strpos(f_base,'.',/reverse_search)

if (o_name eq '') then begin
  if((n_dot le 0) or (n_base-n_dot ne 4)) then begin
    o_name=f_base+'_sieved.log'
  endif else begin
    o_name=strmid(f_base,0,n_dot)+'_sieved.log'
  endelse
endif

out_name=dialog_pickfile(title='Select output point cloud file name',file=o_name, $
dialog_parent=event.top,path=f_path)

;now check for no name, output=input or existing file
if(out_name eq '') then begin
  result=dialog_message('Try again ? (No/Yes) ',/question,/default_no, $
  title='No output file selected!',$
  dialog_parent=event.top)
  if(result eq 'No') then begin
    goto, cleanup
  endif else goto, output
;We will not let the output be the input file
endif else if (strtrim(out_name,2) eq $
  strtrim(cal_name,2)) then begin
  info_text=[strtrim('File name conflict',2),$
             strtrim('Output cannot be the Input !',2),$
             strtrim('Try another selection',2)]
  result=dialog_message(info_text,/error,title='Invalid Input',$
  dialog_parent=event.top)
  out_name=''
  o_name=''
  goto, output
endif

;see if the output file exists
if(file_test(out_name)) then begin
  result=dialog_message('OK to overwrite '+strtrim(out_name,2)+' ? (Yes/No) ', $
  /question,/default_no, $
  title='Output ptcl File Exists',$
  dialog_parent=event.top)
  if(result eq 'Yes') then begin
    fids=envi_get_file_ids()
    if(fids[0] eq -1) then begin
      file_delete, out_name
    endif else begin
      for i=0,n_elements(fids)-1 do begin
        envi_file_query,fids[i],fname=tname
        if (strtrim(out_name,2) eq $
            strtrim(tname,2)) then begin
            envi_file_mng,id=fids[i],/remove
        endif
      endfor
      file_delete, out_name
    endelse
  endif else begin
    out_name=''
    o_name=''
    goto, output
  endelse
endif

;print,'out_name=',strtrim(out_name,2)

text_err=0
openw, tfile, out_name,/get_lun,error=text_err
if (text_err ne 0) then begin
  info_text=[strtrim(string('Error opening log file !'),2),$
    strtrim('File Name ='+strtrim(out_name,2),2),$
    strtrim('Error Number ='+strtrim(text_err,2),2),$
    strtrim('Error Type ='+strtrim(string(!ERROR_STATE.MSG),2),2)]
  result=dialog_message(info_text,/error,title='Error opening log file in DWEL_point_cloud')
  print,'Error opening log file in point cloud!!'
  print,'File Name =',strtrim(out_name,2)
  print,'text_err=',text_err
  print,'Error Type =',strtrim(string(!ERROR_STATE.MSG),2)
  print,'Sys_Error Type =',strtrim(string(!ERROR_STATE.SYS_MSG),2)

  error=1b
  goto, cleanup
endif

time_date=strtrim(systime(),2)

printf,tfile,strtrim('[Point Cloud Data]',2)
printf,tfile,strtrim('Run made at: '+time_date,2)
;printf,tfile,strtrim('X,Y,Z,Intensity,Return_Number,Number_of_Returns,range,theta,phi,Sample,Line,Band',2)
printf,tfile,strtrim('X,Y,Z,d_I,Return_Number,Number_of_Returns,Shot_Number,Run_Number,range,theta,phi,Rk,Sample,Line,Band',2)
flush,tfile

widget_control,/hourglass

count=0L
for j=0L,n_points-1L do begin
  if(((x[j] ge min_x) and (x[j] le max_x)) and $
    ((y[j] ge min_y) and (y[j] le max_y)) and $
    ((z[j] ge min_z) and (z[j] le max_z)) and $
    ((intensity[j] ge min_i) and (intensity[j] le max_i))) then begin
      buf=string(x[j],y[j],z[j],intensity[j],fix(temp[4,j]),fix(temp[5,j]),long(temp[6,j]),$
        long(temp[7,j]),temp[8,j],temp[9,j],$
        temp[10,j],temp[11,j],fix(temp[12,j]),fix(temp[13,j]),fix(temp[14,j]),format='(3f14.3,f14.4,2i10,2i14,4f14.3,3i10)')
      buf=strtrim(strcompress(buf),2)
      while (((ii = strpos(buf, ' '))) ne -1) do $
        strput, buf, ',', ii
      printf,tfile,buf
      count=count+1L
  endif
endfor

flush,tfile
widget_control
free_lun,tfile,/force

;get path and file name as separate strings
o_base=file_basename(out_name)
o_path=file_dirname(out_name)

info_text=[ $
'Point Cloud File Sieving Complete',$
'',$
'Output File Name: '+strtrim(o_base,2),$
'Path to file: '+strtrim(o_path,2),$
'',$
'Number of points in input cloud: '+strtrim(string(n_points,format='(i10)'),2),$
'Number of points in output cloud: '+strtrim(string(count,format='(i10)'),2),$
'',$
'Run another sieve (Y/N) ?' $
]

result=dialog_message(info_text, $
  /question,/default_no, $
  title='Sieving Finished - Run another?',$
  dialog_parent=event.top)
if(result eq 'Yes') then begin
  goto,input
endif

;print,'number of point going through sieve=',count

;exit to caller
cleanup:
pout=0b
band_pos=0b
gain=0b
offset=0b
accum=0b
free_lun, tfile,/force
;get state pointer
widget_control,event.top,get_uvalue=pstate
;clean up pointers
widget_control,event.top,/destroy

return

end

;======================================================================

function pt_util_bannertext
compile_opt idl2

;Enter the banner text as a string array

text=[ $
'Module Description:',$
' ',$
'Procedure pt_util allows you to sieve ascii (name.log) file point clouds',$
'in x, y, z and Intensity. This is a simple routine to use and is needed since',$
'ENVI Lidar tools do not sieve in z or intensity and also because a reduced',$
'ascii file can often be put into a spreadsheet or manipulated by simple programs',$
'that do not read and/or write LAS files.',$
'',$
'The routine reads the ascii log file and reports its ranges in the values',$
'It then asks you to specifiy the new ranges and writes the new log file',$
'At this time it does not re-write the pch (metadata) file and is not very smart',$
'in its checking - so treat it carefully :-)',$
'' $
]

return,text

end

pro pt_util_action_go, event
compile_opt idl2

;Event handler for hitting the "GO" button
;simply calls the main action routine

widget_control,event.top,get_uvalue=pstate

(*pstate).go=1

pt_util_doit, event

return

end

pro pt_util_action_exit, event
compile_opt idl2

;Event handler for hitting the "EXIT" button
;cleans up and destroys the widget hierarchy

widget_control,event.top,get_uvalue=pstate
ptr_free, pstate
widget_control,event.top,/destroy

;some protective programming for now
heap_gc,/verbose

return

end

pro pt_util_resize, event
compile_opt idl2

;First setup protective aunty-Catch to watch out for errors
error_status=0
catch, error_status
if (error_status ne 0) then begin
   catch,/cancel
   help,/last,output=out
   info_text=[strtrim(string('Catch trapped an Error !'),2),$
              strtrim('Error Name : '+strtrim(!err_string,2),2),$
              strtrim(string('Error Number: ',error_status,$
              format='(a,i5)'),2),$
              strtrim('Last Message: '+strtrim(out,2),2)]
   result=dialog_message(info_text,/error,title='Error in pt_util_resize', $
   dialog_parent=event.top)
   goto, cleanup
endif

;Window resize event handler

 widget_control,event.top,get_uvalue=pstate
 widget_control,event.top, xsize=event.x,ysize=event.y

 xs=max([53*(event.x/(*pstate).xy_old[0]),1])
 ys1=max([10*(event.y/(*pstate).xy_old[1]),1])

 if((event.y/(*pstate).xy_old[1]) eq 1) then begin
  ys2=1
 endif else begin
  ys2=max([2*(event.y/(*pstate).xy_old[1]),1])
 endelse

 xs_act=max([(*pstate).xy_act[0]*(event.x/(*pstate).xy_old[0]),1])
 ys_act=min([max([(*pstate).xy_act[1]*(event.y/(*pstate).xy_old[1]),1]),35])

 widget_control,(*pstate).w_2,xsize=xs,ysize=ys1
 widget_control,(*pstate).w_3,xsize=xs,ysize=ys2
 widget_control,(*pstate).wb_action,xsize=xs_act,ysize=ys_act

 widget_control,event.top,/realize

cleanup:

return

end

pro pt_util_cleanup,top
compile_opt idl2

;The cleanup routine for Xmanager

;get state pointer
widget_control,top,get_uvalue=pstate

;clean up pointers
ptr_free, pstate

return

end

pro pt_util, event
compile_opt idl2

;The main routine that sets up the widget hierarchy and manages the events
;There are only two main events here so it is not so hard

; Setup aunty-Catch to watch out for nasty errors
error_status=0
catch, error_status
if (error_status ne 0) then begin
   catch,/cancel
   help,/last,output=out
   info_text=[strtrim(string('Catch trapped an Error !'),2),$
              strtrim('Error Name : '+strtrim(!err_string,2),2),$
              strtrim(string('Error Number: ',error_status,$
              format='(a,i5)'),2),$
              strtrim('Last Message: '+strtrim(out,2),2)]
   result=dialog_message(info_text,/error,title='Error in pt_util')
   goto, cleanup
endif

xoff=200
yoff=150

;now set up the top level base widget
envi_center, xoff,yoff
wb_base=widget_base(title='Sieving Point Cloud data',/column,$
xoffset=xoff,yoffset=yoff,/frame,$
/tlb_size_events)

rev=2s

;get the procedure description and disclaimer to display
intro_text=pt_util_bannertext()
dis_claim=DWEL_disclaimer(rev)

;Set up a base widget to hold the description and disclaimer
wb_info=widget_base(wb_base,/column)

;Use two slider label box routines from the ENVI library
w_2=widget_slabel(wb_info,prompt=intro_text,/frame,xsize=53,ysize=10)
w_3=widget_slabel(wb_info,prompt=dis_claim,/frame,xsize=53,ysize=1)

;Now set up another base widget for the action buttons
wb_action=widget_base(wb_base,/row)

;Action buttons are just "go" and "exit"
w_go=widget_button(wb_action,value='Go',uvalue='pt_util_action_go',$
event_pro='pt_util_action_go',/dynamic_resize)
w_exit=widget_button(wb_action,Value='Exit',uvalue='pt_util_action_exit',$
event_pro='pt_util_action_exit',/dynamic_resize)

;realise the widget hierarchy
widget_control,wb_base,/realize

go=0

widget_geometry_1=widget_info(wb_base,/geometry)
widget_geometry_2=widget_info(wb_action,/geometry)

xy_old=[widget_geometry_1.xsize,widget_geometry_1.ysize]
xy_act=[widget_geometry_2.xsize,widget_geometry_2.ysize]

;set up the state structure
state={ $
tlb:wb_base,$
w_2:w_2, $
w_3:w_3, $
wb_action:wb_action, $
go:go, $
xy_old:xy_old, $
xy_act:xy_act, $
wb_info:wb_info $
}

pstate=ptr_new(state,/no_copy)
widget_control,wb_base,set_uvalue=pstate

;call xmanager
xmanager,'pt_util',wb_base, $
event_handler='pt_util_resize',$
cleanup='pt_util_cleanup',/no_block

;exit back to ENVI
;leave a bit of a trace for bad exits for now
cleanup:

heap_gc,/verbose

return

end
