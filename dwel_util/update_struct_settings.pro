function update_struct_settings, defaultsettings, settings
;+
;PURPOSE:
;; Generate a structure of updated parameters in tunable settings of a given
;; DWEL processing program from the given default settings and user-defined
;; settings. 
;;
;INPUTS:
;; defaultsettings = structure, one pair of tag and value gives a default value
;; of one setting, e.g. {xmin:-50.0, xmax:+50.0}
;;
;; settings = structure, one pair of tag and value gives a user-defined value of
;; one setting. The tag name in this structure must be one of the tag names of
;; defaultsettings. Otherwise, default setting value will be used and print out
;; a warning. If no value of a setting is provided, default setting value is
;; used. 
;;
;OUPUTS:
;; None.
;;
;RETURN:
;; structure, the final values of all settings updated by user-defined values. 
;;
;MODIFICATION HISTORY:
;; Created, 20150121, by Zhan Li, zhanli86@bu.edu
;-
;

  finalsettings = defaultsettings
  ;; tag names we need in settings
  setting_tag_names = tag_names(finalsettings)
  if n_elements(settings) ne 0 or arg_present(settings) then begin
    ;; user supplied settings
    ;; go through user supplied settings and update the final settings. 
    numtags = n_tags(settings)
    tags = tag_names(settings)
    for n=0, numtags-1 do begin
      tmpind = where(strmatch(setting_tag_names, tags[n], /fold_case) eq 1, $
        tmpnum) 
      if tmpnum eq 1 then begin
        finalsettings.(tmpind) = settings.(n)
      endif else begin
        print, 'Tag name is invalid or ambiguous in given ' + $
          'settings. The following setting value is not updated. '
        print, 'Given tag name = ' + strtrim(tags[n], 2)
        print, 'Check the tag name of this setting.'
      endelse 
    endfor 
  endif 

  return, finalsettings

end
