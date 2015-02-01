%% Extract panel returns from stationary scans
% Clean point cloud of stationary scans and extract those points actually
% from the panel. Then point cloud is generated from DWEL2.0 program by
% faking a data cube from stationary scans. 
% 
% Zhan Li, zhanli86@bu.edu
% Created: 20141104
% Last modified: 20141104

if ~exist('enviinfo.m', 'file')
    addpath(genpath('~/Documents/MATLAB/envi-rw'), path);
end

% inptclfile = ['/projectnb/echidna/lidar/DWEL_Processing/DWEL_TestCal/cal-nsf-' ...
%               '20140812/7/' ...
%               'cal_nsf_20140812_7_1064_cube_bsfix_pxc_update_ptcl_points.txt'];
% outpanelreturnsfile = ['/projectnb/echidna/lidar/DWEL_Processing/' ...
%                     'DWEL_TestCal/cal-nsf-20140812/7/' ...
%                     'cal_nsf_20140812_7_1064_cube_bsfix_pxc_update_ptcl_points_panel.txt'];
% satinfofile = ['/projectnb/echidna/lidar/DWEL_Processing/DWEL_TestCal/cal-' ...
%                'nsf-20140812/7/' ...
%                'cal_nsf_20140812_7_1064_cube_bsfix_pxc_update_ancillary.img'];
% panelrange = 7;

% read the point cloud
fid = fopen(inptclfile, 'r');
if fid == -1 
    fprintf('Error openning input point cloud file: %s\n', inptclfile)
    return
end
data = textscan(fid, repmat('%f', 1, 16), 'HeaderLines', 3, 'Delimiter', ',');
fclose(fid);
intensity = data{4};
returnind = data{5};
numreturns =  data{6};
shotnum = data{7};
range = data{9};
sampleind = data{13};
lineind = data{14};
p_int = data{4};
p_fwhm = data{16};
clear data;

% read saturation mask
hdrinfo = envihdrread([satinfofile, '.hdr']);
satmask = multibandread(satinfofile, [hdrinfo.lines, hdrinfo.samples, ...
                    hdrinfo.bands], hdrinfo.data_type.matlab, hdrinfo.header_offset, ...
                        hdrinfo.interleave, hdrinfo.byte_order.matlab, {'Band', ...
                   'Direct', '1'});

% remove some outliers from side lobes or noise
tmpind = range > panelrange-1 & range < panelrange+1 & returnind == 1 & numreturns <= 2;
% tmpind = range > 6.95 & range < 7.1 & returnind == 1 & numreturns <= 2;
% get saturation mask value for each point. 
linearind = sub2ind(size(satmask), lineind, sampleind);
satmaskval= satmask(linearind);

% intensity = intensity(tmpind);
% range = range(tmpind);
% sampleind = sampleind(tmpind);
% lineind = lineind(tmpind);
panelreturns = [shotnum(tmpind), intensity(tmpind), numreturns(tmpind), range(tmpind), ...
                p_int(tmpind), p_fwhm(tmpind), sampleind(tmpind), lineind(tmpind), satmaskval(tmpind)];

% sort out the points, ascending lineind first and ascending sampleind second
panelreturns = sortrows(panelreturns, [8, 7]);

% write out
fid = fopen(outpanelreturnsfile, 'w');
fprintf(fid, 'shot_num,d_out,number_of_returns,range,I,FWHM,sample,line,sat_mask\n');
fprintf(fid, ['%d,%.3f,%d,%.3f,%.3f,%.3f,%d,%d,%d\n'], panelreturns');
fclose(fid); 
