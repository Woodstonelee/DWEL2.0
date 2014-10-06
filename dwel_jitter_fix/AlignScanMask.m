function AlignedWfMax_dual = AlignScanMask(datacubefile1064, datacubefile1548, gapopen, startline, AlignedMaskFile)
    % parameters initialization
    elements={'samples' 'lines' 'bands' 'data type'};
    d={'bit8' 'int16' 'int32' 'float32' 'float64' 'uint16' 'uint32' 'int64' 'uint64'};

    % Open ENVI header file to retreive s, l, b & d variables
    rfid = fopen(strcat(datacubefile1064,'.hdr'),'r');

    % Check if the header file is correctely open
    if rfid == -1
        error('Input header file does not exist');
    end

    % Read ENVI image header file 
    while 1
        tline = fgetl(rfid);
        if ~ischar(tline), break, end
        [first,second]=strtok(tline,'=');
        first = strtrim(first);
        second = strtok(second, '=');
        second = strtrim(second);
        
        switch first
          case elements(1)
            ns=str2num(second);
          case elements(2)
            nl=str2num(second);
          case elements(3)
            nb=str2num(second);
          case elements(4)
            t=str2num(second);
            switch t
              case 1
                t=d{1};
              case 2
                t=d{2};
              case 3
                t=d{3};
              case 4
                t=d{4};
              case 5
                t=d{5};
              case 12
                t=d{6};
              case 13
                t=d{7};
              case 14
                t=d{8};
              case 15
                t=d{9};
              otherwise
                error('Unknown image data type');
            end
        end
    end
    fclose(rfid);
        
    WfMax_1064 = multibandread(datacubefile1064, [nl, ns, nb], t, 0, 'bil', 'ieee-le', {'Band', 'Direct', 6});
    WfMax_1548 = multibandread(datacubefile1548, [nl, ns, nb], t, 0, 'bil', 'ieee-le', {'Band', 'Direct', 6});

    WfMax_dual = zeros(nl, ns, 2);
    WfMax_dual(:,:,1) = WfMax_1064(1:nl,:);
    WfMax_dual(:,:,2) = WfMax_1548(1:nl,:);
    
    % nc = 128;
    % tmpmin = min(min(WfMax_dual(:,:,1)));
    % tmp = (WfMax_dual(:,:,1)-tmpmin)/(max(max(WfMax_dual(:,:,1))) ...
    %                                          - tmpmin)*nc;
    % imwrite(tmp, jet(nc), 'test.tif', 'tif');
    % AlignedWfMax_dual=WfMax_dual;
    % fprintf('test finish\n');
    
    % gapopen = -1024^2*2;
    [AlignedWfMax_dual, AlignedMask] = NWAlignScanImage_nd(WfMax_dual, gapopen, startline);
    k = strfind(AlignedMaskFile, '.tif');    
    if isempty(k)
        tmp1064file = [AlignedMaskFile, '-', num2str(gapopen), ...
                       '-1064.tif'];
        tmp1548file = [AlignedMaskFile, '-', num2str(gapopen), ...
                       '-1548.tif'];
        AlignedMaskFile = [AlignedMaskFile, '-', num2str(gapopen), '.tif'];
    else 
        tmp1064file = [AlignedMaskFile(1:k(end)-1), '-', num2str(gapopen), ...
                       '-1064.tif'];
        tmp1548file = [AlignedMaskFile(1:k(end)-1), '-', num2str(gapopen), ...
                       '-1548.tif'];
        AlignedMaskFile = [AlignedMaskFile(1:k(end)-1), '-', num2str(gapopen), ...
                       '.tif'];
    end

    imwrite(AlignedMask, AlignedMaskFile, 'tif');
    
    nc = 128;
    tmpmin = min(min(AlignedWfMax_dual(:,:,1)));
    tmp = (AlignedWfMax_dual(:,:,1)-tmpmin)/(max(max(AlignedWfMax_dual(:,:,1))) ...
                                             - tmpmin)*nc;
    imwrite(tmp, jet(nc), tmp1064file, 'tif');
    tmpmin = min(min(AlignedWfMax_dual(:,:,2)));
    tmp = (AlignedWfMax_dual(:,:,2)-tmpmin)/(max(max(AlignedWfMax_dual(:,:,2))) ...
                                             - tmpmin)*nc;
    imwrite(tmp, jet(nc), tmp1548file, 'tif');
end