;+
; NAME:
;DWEL_STATIC_HDF2CUBE
;
;
; PURPOSE:
;Fake a data cube from a stationary scan for DWEL program to process easily
;later. 
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;
;
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;-

function CheckStatDWEL, DWEL_H5File, Wavelength, nsample, err
  compile_opt idl2
  
  scanenc_ind = 0
  rotateenc_ind = 1
  err=0
  
  ;; open the HDF5 file
  fileid=h5f_open(DWEL_H5File)
  
  Waveset_Name = '/'+strtrim(string(Wavelength), 2)+' Waveform Data'
  waveset = h5d_open(fileid, Waveset_Name)
  
  wave_type = h5d_get_type(waveset)
  wave_class = h5t_get_class(wave_type)
  if strcmp(wave_class, 'H5T_INTEGER', /fold_case) then begin
    wavespace = h5d_get_space(waveset) ; space is the actual place where the data can be read out.
    tmpdims = H5S_GET_SIMPLE_EXTENT_DIMS(wavespace)
    NoSamplesPerShot = tmpdims[0] ; number of bins per each waveform.
    numwfs = tmpdims[1] ; number waveforms in this stationary scan
  endif else begin
    print, Waveset_Name + ': data type is not H5T_INTEGER as expected. Please check the HDF5 file.'
    h5t_close, wave_type
    h5d_close, waveset
    h5f_close, fileid
    H5_CLOSE
    err=2
    return, -1
  endelse
  h5t_close, wave_type
  
  h5d_close, waveset
  h5f_close, fileid
  H5_CLOSE  ; This is the crucial step - release all of HDF5's memory

  ;; fabricate all ancillary information for a fake data cube
  nadirelevshift = 262144L
  TotalNoScans = fix(numwfs/nsample) + 1
  NoShotsPerScan = nsample
  shotstart = 0
  shotend = numwfs-1
  ShotStartVec = indgen(TotalNoScans, /long)*NoShotsPerScan
  ShotEndVec = (indgen(TotalNoScans, /long)+1)*NoShotsPerScan - 1
  ShotEndVec[TotalNoScans-1] = numwfs - 1
  ShotNumVec = ShotEndVec - ShotStartVec + 1
  NoScanPerRotation = round(2*!pi/2e-3)

  return, {DWELFileName:DWEL_H5File, NadirScanEncoder:nadirelevshift, $
    TotalNoScans:TotalNoScans, NoShotsPerScan:NoShotsPerScan, NoSamplesPerShot:NoSamplesPerShot, $
    FirstShotInd:shotstart, LastShotInd:shotend, $
    ShotStart:ShotStartVec, ShotEnd:ShotEndVec, ShotNum:ShotNumVec, $
    NoScanPerRotation:NoScanPerRotation}
end

function FakeDataCube, DWEL_MetaInfo, DataCube_File, Wavelength, AncillaryFile
  compile_opt idl2
  
  scanenc_ind = 0
  rotateenc_ind = 1
  AncillaryFile=''
  ;; create the name of the ancillary file from the given waveform
  ;;data cube file name.
  AncillaryFile = strtrim(strmid(DataCube_File,0,strpos(DataCube_File, '.', /reverse_search))+'_ancillary.img',2)
  
  fileid=h5f_open(DWEL_MetaInfo.DWELFileName)
  
  ;; open the dataset of waveform and get the space of the dataset.
  Waveset_Name = '/'+strtrim(string(Wavelength), 2)+' Waveform Data'
  waveset = h5d_open(fileid,Waveset_Name)
  wave_space = h5d_get_space(waveset)
  memspace = h5s_create_simple([DWEL_MetaInfo.NoSamplesPerShot, 1])
  
  openw, DataCubeFID, DataCube_File, /get_lun
  openw, AncillaryFID, AncillaryFile, /get_lun
  
  DataArray = intarr(DWEL_MetaInfo.NoShotsPerScan, DWEL_MetaInfo.NoSamplesPerShot)
  AncillaryArray = lonarr(DWEL_MetaInfo.NoShotsPerScan,9)
  Trigger = 0B
  SunSensor = 0B
  ScanEncoder = 0L
  RotaryEncoder = 0L
  LaserPower = 0E
  Pos = 0L
  mask=0b
  ShotZen=0.0
  ShotAzim=0.0
  
  shotind = DWEL_MetaInfo.FirstShotInd
  NumBlank = 0
  BlankVec = bytarr(DWEL_MetaInfo.NoShotsPerScan)
  ZeroWaveform = intarr(DWEL_MetaInfo.NoSamplesPerShot)
  for i = 0, DWEL_MetaInfo.TotalNoScans-1, 1 do begin
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; if the number of shots in a scan line is less than the given
    ;; NoShotsPerScan (number of column in a row), the following
    ;;generates a vector BlankPixelLoc which gives where we need some
    ;;dummy pixels to fill the extra space in this row.
    ;; check how many blank pixels there will be. Insert the blank
    ;; pixels to fill the gap.
    NumBlank = DWEL_MetaInfo.NoShotsPerScan - (DWEL_MetaInfo.ShotEnd[i] - $
      DWEL_MetaInfo.ShotStart[i] + 1)
    if NumBlank gt 0 then begin
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; simply insert blank pixels at the end of each scan line
      BlankPixelLoc = indgen(NumBlank) + DWEL_MetaInfo.NoShotsPerScan-NumBlank
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    endif else begin ; if no blank pixel is needed
      BlankPixelLoc=-1
    endelse
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    for j = 0L, DWEL_MetaInfo.NoShotsPerScan-1, 1  do begin
      if (total((BlankPixelLoc eq j)) gt 0) or (shotind ge DWEL_MetaInfo.LastShotInd+1) then begin
        mask = 0L
        AncillaryArray[j,*] = [fix(Trigger,type=3), fix(SunSensor,type=3), 0L, 0L, $
          0L, 0L, fix(mask,type=3), 0L, $
          0L]
        DataArray[j,*] = ZeroWaveform
        continue
      endif
      H5S_SELECT_HYPERSLAB, wave_space, [[shotind], [0]], [[1], $
        [DWEL_MetaInfo.NoSamplesPerShot]], /reset
      Waveform = h5d_read(waveset, file_space=wave_space, MEMORY_SPACE = $
        memspace)
      ;; As of 20140913, one waveform from NSF DWEL has first two samples storing
      ;; shot sequence number and ??? number according to Kuravi. Simply
      ;; overwrite them here. May need to put these two numbers into ancillary
      ;; file if they are found to be informative later.
      Waveform[0] = Waveform[2]
      Waveform[1] = Waveform[2]
      
      WaveformMean = mean(waveform)
      WaveformMax = max(waveform, max_I)
      DataArray[j,*] = Waveform
      ;; fake encoder values
      ScanEncoder = fix(524288*(1 - $
        (j+0.5)/double(DWEL_MetaInfo.NoShotsPerScan)), type=3)
      RotaryEncoder = fix(524288*(1 - (i+0.5)*2e-3/(2*!pi)), type=3)
      
      if ((RotaryEncoder eq 0) or (ScanEncoder eq 0)) then begin
        RotaryEncoder=0L
        ScanEncoder=0L
        ShotZen=0.0
        ShotAzim=0.0
        mask=0L
        WaveformMean = 0L
        WaveformMax = 0L
        DataArray[j,*] = ZeroWaveform
      endif else begin
        ShotZen = ((double(DWEL_MetaInfo.NadirScanEncoder) - double(ScanEncoder)) / double(524288)) * 360.0d0
        ShotAzim = ((double(524288)-double(RotaryEncoder))/double(524288)) * 360.0d0
        
        if (ShotZen lt -180.0) then ShotZen=ShotZen+360.0
        if (ShotZen gt 180.0) then ShotZen = ShotZen - 360.0
        
        if (ShotZen lt 0.0) then begin
          ShotZen=-ShotZen
          ShotAzim=ShotAzim+180.0
        endif
        if (ShotAzim gt 360.0) then ShotAzim=ShotAzim-360.0
        
        mask = 1L
        
      endelse
      AncillaryArray[j,*] = [fix(Trigger,type=3), fix(SunSensor,type=3), ScanEncoder, RotaryEncoder, $
        fix(round(LaserPower*100.0),type=3), fix(round(WaveformMax),type=3), fix(mask,type=3), fix(round(100.0*ShotZen),type=3), $
        fix(round(100.0*ShotAzim),type=3)]
      shotind = shotind + 1L
    endfor
    
    writeu, DataCubeFID, DataArray
    writeu, AncillaryFID, AncillaryArray
  endfor
  
  h5s_close, wave_space
  h5s_close, memspace
  h5d_close, waveset
  h5f_close, fileid
  ; This is the crucial step - release all of HDF5's memory
  H5_CLOSE
  
  close, DataCubeFID
  close, AncillaryFID
  free_lun, DataCubeFID, /force
  free_lun, AncillaryFID, /force
  
  data_dims = size(DataArray)
  ancillary_dims = size(AncillaryArray)
  
  return, {samples:data_dims[1], lines:DWEL_MetaInfo.TotalNoScans, $
    databands:data_dims[2], ancillarybands:ancillary_dims[2], offset:0, $
    filetype:'ENVI Data Cube', datatype:data_dims[3], $
    ancillarydatatype:ancillary_dims[3], interleave:1, sensortype:'DWEL', $
    byteorder:0, wavelengthunit:'metres', range:120.0}
    
end

pro dwel_static_hdf2cube, instathdf5file, Config_File, outcubefile, wavelength, $
  wavelength_label, DWEL_Height, beam_div, srate, nsample=nsample

  compile_opt idl2

  ;; check the keyword argument
  if n_elements(nsample) ne 0 or arg_present(nsample) then begin
    ;; user gives a nsample
  endif else begin
    ;; no nsample given, use default values
    nsample = 3142
  endelse 

  inlun=30
  anc_name=''
  out_fid=0
  anc_fid=0
  AncillaryFile=''
  err_flag=0
  err=0

  ;; get the size of input file to be processed. It will be used in later
  ;; summary of processing time. 
  procfilesize = file_info(instathdf5file)
  procfilesize = procfilesize.size
  ;; get the time now as the start of processing
  starttime = systime(1)
    
  DWEL_MetaInfo = CheckStatDWEL(instathdf5file,Wavelength,nsample,err)
  
  if (err ne 0) then begin
    print,'error return from CheckStatDWEL'
    print,'Local code='+strtrim(string(err),2)
    err_flag=1
    goto,out
  endif
  HeaderInfo = FakeDataCube(DWEL_MetaInfo, outcubefile, Wavelength, AncillaryFile)
  
  ;get path and DWEL_file name as separate strings
  last=strpos(instathdf5file,path_sep(),/reverse_search)
  f_path=strmid(instathdf5file,0,last+1)
  f_base=strtrim(strmid(instathdf5file,last+1,strlen(instathdf5file)-last-1),2)
  seek=strlowcase(strtrim(f_base,2))
  i1=strpos(seek,'_')
  i2=strpos(seek,'.hdf5')
  date_time=strmid(seek,i1+1,i2-i1-1)
  print,'Date_Time=',date_time

  istat=dwel_get_config_info(Config_File,Config_Info,consum)
  if (istat gt 0) then begin
    print,'error in call to dwel_get_config_info'
    consum=''
    config_info=['']
    err_flag=2
    goto,out
  endif
  print,'consum=',consum
    
  Name_Info=['Program=dwel2cube_cmd_nsf',$
    'Original DWEL HDF File='+strtrim(f_base,2)]
  Site_Info=[ $
    'DWEL Date Time='+strtrim(date_time,2),$
    'Processing Date Time='+strtrim(systime(),2)]
  Scan_Info=[ $
    'Beam Divergence='+strtrim(string(beam_div,format='(f10.2)'),2)+' mrad',$
    'Number of Hemispheres='+strtrim(string(1),2),$
    'Scans per Complete Rotation='+strtrim(string(DWEL_MetaInfo.NoScanPerRotation),2),$
    'Number of Scans selected='+strtrim(string(DWEL_MetaInfo.TotalNoScans),2),$
    'Number of Shots per Scan='+strtrim(string(DWEL_MetaInfo.NoShotsPerScan),2),$
    'Number of Samples per Shot='+strtrim(string(DWEL_MetaInfo.NoSamplesPerShot),2),$
    'Digitised Range='+strtrim(string(120.0,format='(f10.2)'),2)+' Metres',$
    'Digitiser Sampling Rate='+strtrim(string(srate,format='(f10.2)'),2)+' smp/ns']
  Post_Info=[ $
    'Data Start='+strtrim(string(0),2),$
    'Actual scans completed='+strtrim(string(DWEL_MetaInfo.TotalNoScans),2)]
    
  DWEL_Scan_Info=[Name_Info,Site_Info,Scan_Info,Post_Info,Config_Info]
  
  if (DWEL_Height ge 0) then $
    DWEL_Scan_Info=[DWEL_Scan_Info,'DWEL Height='+strtrim(DWEL_Height,2)]
    
  DSR=srate
  wl_out=findgen(HeaderInfo.databands)/DSR
  band_names=strarr(HeaderInfo.databands)+'Time_Sample_'
  band_names=band_names+strtrim(string(indgen(HeaderInfo.databands)+1),2)
  
  ;get output_file name without path
  last=strpos(outcubefile,path_sep(),/reverse_search)
  out_base=strtrim(strmid(outcubefile,last+1,strlen(outcubefile)-last-1),2)
  DWEL_Scan_Info=[DWEL_Scan_Info,'Output Cube File='+out_base]
  
  DWEL_Adaptation=['Band "Waveform Mean" is actually "Waveform Max"', 'Band "Scan Encoder" is value corrected for nadir shift']
  DWEL_Adaptation=[DWEL_Adaptation, 'Wavelength='+strtrim(Wavelength_Label, 2)]
  DWEL_Adaptation=[DWEL_Adaptation, 'Scan encoder of zenith point='+strtrim(DWEL_MetaInfo.NadirScanEncoder, 2)]
  
  ENVI_SETUP_HEAD, fname=outcubefile, $
    ns=HeaderInfo.samples, nl=HeaderInfo.lines, nb=HeaderInfo.databands, $
    interleave=HeaderInfo.interleave, data_type=HeaderInfo.datatype, $
    offset=HeaderInfo.offset, zplot_titles=['Time (nsec)','Intensity'], $
    bnames=band_names, $
    wl=wl_out, /write
    
  envi_open_file, outcubefile,r_fid=out_fid,/no_interactive_query,/no_realize
  
  tname=''
  envi_file_query,out_fid,fname=tname
  print,'cube file name='+strtrim(tname,2)
  
  envi_assign_header_value, fid=out_fid, keyword='DWEL_Scan_Info', $
    value=DWEL_Scan_Info
  envi_assign_header_value, fid=out_fid, $
    keyword='DWEL_Adaptation', $
    value=DWEL_Adaptation
  envi_write_file_header, out_fid
  
  envi_file_mng, id=out_fid,/remove
  
  bnames=['Non Triggers','Sun Sensor','Scan Encoder','Rotary Encoder', $
    'Laser Power','Waveform Mean','Mask','Zenith','Azimuth']
  ENVI_SETUP_HEAD, fname=AncillaryFile, $
    ns=HeaderInfo.samples, nl=HeaderInfo.lines, nb=HeaderInfo.ancillarybands, $
    interleave=HeaderInfo.interleave, data_type=HeaderInfo.ancillarydatatype, $
    offset=HeaderInfo.offset,bnames=bnames,/write
    
  envi_open_file, AncillaryFile,r_fid=anc_fid,/no_interactive_query,/no_realize
  
  tname=''
  envi_file_query,anc_fid,fname=tname
  print,'anc file name='+strtrim(tname,2)
  
  envi_assign_header_value, fid=anc_fid, $
    keyword='DWEL_Scan_Info', $
    value=DWEL_Scan_Info
  envi_assign_header_value, fid=anc_fid, $
    keyword='DWEL_Adaptation', $
    value=DWEL_Adaptation
    
  envi_write_file_header, anc_fid
  
  envi_file_mng, id=anc_fid,/remove
  
  out:
  envi_file_mng, id=out_fid,/remove
  envi_file_mng, id=anc_fid,/remove
  free_lun,inlun,/force
  if (err_flag ne 0) then print,'Returning with error from dwel2cube_cmd_nsf'
  heap_gc,/verbose
  
  ;; write processing time summary
  print, '******************************************'
  print, 'Processing program = dwel_stationary2cube'
  print, 'Input HDF5 file size = ' + $
    strtrim(string(double(procfilesize)/(1024.0*1024.0*1024.0)), 2) + ' G'
  print, 'Processing time = ' + strtrim(string((systime(1) - starttime)), 2) + ' ' + $
    'seconds'
  print, '******************************************'

  return
  
end
