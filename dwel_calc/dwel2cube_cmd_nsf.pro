; Convert DWEL HDF5 data format to ENVI data cube (.img files).
; Running in command line mode.
;
; Zhan Li, zhanli86@bu.edu
; Adapted from heritage EVI code, 2012
; Revison: March 2014

function CheckDWEL, DWEL_H5File, Wavelength, nadirelevshift, err
  compile_opt idl2
  
  scanenc_ind = 0
  rotateenc_ind = 1
  err=0
  
  ;; open the HDF5 file
  fileid=h5f_open(DWEL_H5File)
  ;; read encoder dataset.
  encoderset = h5d_open(fileid,'/Interpolated angles (Alt, Azm)')
  encoders = h5d_read(encoderset)
  dim_encoders = size(encoders, /dimensions)
  
  ;; use the difference between every two scan encoder values to
  ;; determine the actual start and ending shots and to remove the
  ;;dummy shots at the beginning and the end.
  interval_diff = long(encoders[scanenc_ind, 0:dim_encoders[1]-2] - encoders[scanenc_ind, 1:dim_encoders[1]-1]) ; the difference between two consecutive shots, the early one - the later one
  tmpind = where(interval_diff ne 0, tmpcount, ncomplement=count, complement=dummyind)
  if (tmpcount gt 0) then begin
    shotstart = long(tmpind[0])
    shotend = long(tmpind[size(tmpind, /n_elements)-1])
  endif else begin
    print, 'No valid scan encoder value! Processing is terminated!'
    h5d_close, encoderset
    h5f_close, fileid
    h5_close
    err=1
    return, -1
  endelse
  
  ;; Find out the number of shots in each scan by identifying the
  ;;start and ending of a scan line.
  ;; If the difference is negative (actually comparing the difference
  ;; with -2^18 instead of zero to avoid possible positive differences
  ;;due to wiggles in the encoder values), it indicates the start of a
  ;;new scan.
  NegInd = where(interval_diff lt -262144, NegCount)
  TotalNoScans = NegCount + 1L   ; TotalNoScans is the dimension of the azimuth axis
  ShotStartVec = [shotstart, NegInd+1L] ; the indexes of start positions of each scan line
  ShotEndVec = [NegInd, shotend] ; the indexes of ending positions of each scan line
  ShotNumVec = ShotEndVec - ShotStartVec + 1L ; the number of shots per each scan line
  
  ;; number of scan lines per a whole 360-degree rotation.
  NoScanPerRotation = long(TotalNoScans * 524288.0 / abs(encoders[rotateenc_ind,0]-encoders[rotateenc_ind,shotend]))
  
  ;; Get the largest number of shots in a scan. It is the dimension of the zenith axis
  NoShotsPerScan = max(ShotNumVec)
  
  Waveset_Name = '/'+strtrim(string(Wavelength), 2)+' Waveform Data'
  waveset = h5d_open(fileid, Waveset_Name)
  
  wave_type = h5d_get_type(waveset)
  wave_class = h5t_get_class(wave_type)
  if strcmp(wave_class, 'H5T_INTEGER', /fold_case) then begin
    wavespace = h5d_get_space(waveset) ; space is the actual place where the data can be read out.
    tmpdims = H5S_GET_SIMPLE_EXTENT_DIMS(wavespace)
    NoSamplesPerShot = tmpdims[0] ; number of bins per each waveform.
  endif else begin
    print, Waveset_Name + ': data type is not H5T_INTEGER as expected. Please check the HDF5 file.'
    h5t_close, wave_type
    h5d_close, waveset
    h5d_close, encoderset
    h5f_close, fileid
    H5_CLOSE
    err=2
    return, -1
  endelse
  h5t_close, wave_type
  
  h5d_close, waveset
  h5d_close, encoderset
  h5f_close, fileid
  H5_CLOSE  ; This is the crucial step - release all of HDF5's memory
  return, {DWELFileName:DWEL_H5File, NadirScanEncoder:nadirelevshift, $
    TotalNoScans:TotalNoScans, NoShotsPerScan:NoShotsPerScan, NoSamplesPerShot:NoSamplesPerShot, $
    FirstShotInd:shotstart, LastShotInd:shotend, $
    ShotStart:ShotStartVec, ShotEnd:ShotEndVec, ShotNum:ShotNumVec, $
    NoScanPerRotation:NoScanPerRotation}
end

function DataCube, DWEL_MetaInfo, DataCube_File, Wavelength, AncillaryFile
  compile_opt idl2
  
  scanenc_ind = 0
  rotateenc_ind = 1
  AncillaryFile=''
  ;; create the name of the ancillary file from the given waveform
  ;;data cube file name.
  AncillaryFile = strtrim(strmid(DataCube_File,0,strpos(DataCube_File, '.', /reverse_search))+'_ancillary.img',2)
  
  fileid=h5f_open(DWEL_MetaInfo.DWELFileName)
  encoderset = h5d_open(fileid, '/Interpolated angles (Alt, Azm)')
  encoders = h5d_read(encoderset)
  
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
    ScanEncoderVec =  long(encoders[scanenc_ind, DWEL_MetaInfo.ShotStart[i]:DWEL_MetaInfo.ShotEnd[i]])
    ;; check how many blank pixels there will be. Insert the blank
    ;; pixels to fill the gap.
    NumBlank = DWEL_MetaInfo.NoShotsPerScan - size(ScanEncoderVec,/n_elements)
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
      ScanEncoder = long(encoders[scanenc_ind, shotind])
      RotaryEncoder = Long(encoders[rotateenc_ind, shotind])
      
      if ((j eq 0) or (RotaryEncoder eq 0) or (ScanEncoder eq 0)) then begin
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
        ;      ;;;;;; temporarily fake azimuth angle b/c of the wrong azimuth encoders.
        ;      ShotAzim = 2*(i+0.5)*0.001/!pi*180.0
        ;      RotaryEncoder = fix(ShotAzim/360.0*double(524288), type=size(ScanEncoder, /type))
        
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
  h5d_close, encoderset
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

pro dwel2cube_cmd_nsf, DWEL_H5File, Config_File, DataCube_File, Wavelength, $
    Wavelength_Label, DWEL_Height, beam_div, srate, err_flag, nadirshift=nadirelevshift
  ;;
  ;; Because the names of the waveform datasets in HDF5 files were
  ;;incorrectly labeled as of March 2014, including all data from CA
  ;;Sierra June 2013 and Brisbane August 2013, two numbers of wavelength
  ;;are required here.
  ;; Wavelength: a number used to create a dataset name and get waveform
  ;; data from HDF5 file.
  ;; Wavelength_Label: a number to be put in the header file. This is
  ;;the CORRECT wavelength number.
    
  compile_opt idl2
  envi, /restore_base_save_files
  envi_batch_init, /no_status_window
  
  inlun=30
  anc_name=''
  out_fid=0
  anc_fid=0
  AncillaryFile=''
  err_flag=0
  err=0
  
  ;; give a default scan encoder value of nadir shift here if you get one so
  ;; that you don't have to put nadir shift in the command everytime
  if n_elements(nadirelevshift) ne 0 or arg_present(k) then begin ; function
  ; calling gives nadir shift value
  
  endif else begin ; nadir shift value is not given, use default values
    nadirelevshift = 392433 ; 130289 ; this default value was from NSF DWEL 2014/05/03
  endelse
  
  DWEL_MetaInfo = CheckDWEL(DWEL_H5File,Wavelength,nadirelevshift,err)
  
  if (err ne 0) then begin
    print,'error return from CheckDWEL'
    print,'Local code='+strtrim(string(err),2)
    err_flag=1
    goto,out
  endif
  HeaderInfo = DataCube(DWEL_MetaInfo, DataCube_File, Wavelength, AncillaryFile)
  
  ;get path and DWEL_file name as separate strings
  last=strpos(DWEL_H5File,path_sep(),/reverse_search)
  f_path=strmid(DWEL_H5File,0,last+1)
  f_base=strtrim(strmid(DWEL_H5File,last+1,strlen(DWEL_H5File)-last-1),2)
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
  
  Name_Info=['Program=dwel2cube_cmd_nsf, import DWEL HDF5 file to ENVI cube image',$
    'Original DWEL HDF File='+strtrim(f_base,2)]
  Site_Info=[ $
    'Scan Description='+strtrim('DWEL_nsf Scan: '+strtrim(consum,2),2),$
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
  last=strpos(DataCube_File,path_sep(),/reverse_search)
  out_base=strtrim(strmid(DataCube_File,last+1,strlen(DataCube_File)-last-1),2)
  DWEL_Scan_Info=[DWEL_Scan_Info,'Output Cube File='+out_base]
  
  DWEL_Adaptation=['Band "Waveform Mean" is actually "Waveform Max"', 'Band "Scan Encoder" is value corrected for nadir shift']
  DWEL_Adaptation=[DWEL_Adaptation, 'Wavelength='+strtrim(Wavelength_Label, 2)]
  DWEL_Adaptation=[DWEL_Adaptation, 'Nadir shift of scan encoder='+strtrim(nadirelevshift, 2)]
  
  ENVI_SETUP_HEAD, fname=DataCube_File, $
    ns=HeaderInfo.samples, nl=HeaderInfo.lines, nb=HeaderInfo.databands, $
    interleave=HeaderInfo.interleave, data_type=HeaderInfo.datatype, $
    offset=HeaderInfo.offset, zplot_titles=['Time (nsec)','Intensity'], $
    bnames=band_names, $
    wl=wl_out, /write
    
  envi_open_file, DataCube_File,r_fid=out_fid,/no_interactive_query,/no_realize
  
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
  if (err_flag ne 0) then print,'Returning with error from dwel2cube_cmd'
  heap_gc,/verbose
end
