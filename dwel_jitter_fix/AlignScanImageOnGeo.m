% Align scan images on geo server

% datacubefile1064 = '/projectnb/echidna/lidar/DWEL_Processing/Brisbane2013Aug/Brisbane2013_Kara05/Aug1_Kara05_C/Aug1_Kara5_C_1064_Cube_ancillary.img';
% datacubefile1548 = '/projectnb/echidna/lidar/DWEL_Processing/Brisbane2013Aug/Brisbane2013_Kara05/Aug1_Kara05_C/Aug1_Kara5_C_1548_Cube_ancillary.img';
% ns = 2703;
% nl = 1670;
% gapopen = -512^2*2;
% startline = 2;
% AlignedMaskFile = '/projectnb/echidna/lidar/DWEL_Processing/Brisbane2013Aug/Brisbane2013_Kara05/Aug1_Kara05_C/Aug1_Kara5_C_AlignMask.tif';

% datacubefile1064 = '/projectnb/echidna/lidar/DWEL_Processing/CA2013June/CA2013_Site305/June14_01_305_NE/June14_01_305_NE_1064_Cube_NadirCorrect_ancillary.img';
% datacubefile1548 = '/projectnb/echidna/lidar/DWEL_Processing/CA2013June/CA2013_Site305/June14_01_305_NE/June14_01_305_NE_1548_Cube_NadirCorrect_ancillary.img';
% ns = 2560;
% nl = 1627;
% gapopen = -1024^2*2;
% startline = 3;
% AlignedMaskFile = '/projectnb/echidna/lidar/DWEL_Processing/CA2013June/CA2013_Site305/June14_01_305_NE/June14_01_305_NE_1064_Cube_NadirCorrect_ancillary_wfmax_AlignedMask.tif';

%fprintf('this is a test\n')
AlignedWfMax_dual = AlignScanMask(datacubefile1064, datacubefile1548, gapopen, startline, AlignedMaskFile);
