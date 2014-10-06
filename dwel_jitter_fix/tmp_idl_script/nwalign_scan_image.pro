pro nwalign_scan_image, anc1064file, anc1548file, alignedmaskfile
  clean_envi_file_fids
  ;read 1064 waveform max intensity image
  envi_open_file, anc1064file, /no_interactive_query, /no_realize, r_fid=anc1064_fid
  if (anc1064_fid eq -1) then begin
    print,strtrim('Error opening input file',2)
    print,'Input File: '+strtrim(anc1064file,2)
    goto, cleanup
  endif
  envi_file_query, anc1064_fid, ns=anc_ns, nl=anc_nl, nb=anc_nb, data_type=anc_type, $
    dims=anc_dims
  ;initialize a 3D array to contain dual-wavelength waveform intensity image
  wfmax_dual = make_array(anc_ns, anc_nl, 2, type=anc_type)
  wfmax_dual[*,*,0] = envi_get_data(dims=anc_dims, fid=anc1064_fid, pos=5)
  ;read 1548 waveform max intensity image
  envi_open_file, anc1548file, /no_interactive_query, /no_realize, r_fid=anc1548_fid
  if (anc1548_fid eq -1) then begin
    print,strtrim('Error opening input file',2)
    print,'Input File: '+strtrim(anc1064file,2)
    goto, cleanup
  endif
  envi_file_query, anc1548_fid, ns=anc_ns, nl=anc_nl, nb=anc_nb, data_type=anc_type, $
    dims=anc_dims
  wfmax_dual[*,*,1] = envi_get_data(dims=anc_dims, fid=anc1548_fid, pos=5)
  
  nd = 2
  gapopen = -512l^2*2
  ;; use new line2 as the next line1 and do nothing to combine two versions of
  ;% alignments.
  ;gfImage = -1*ones(nl, ns*2, nd);
  gfImage = make_array(anc_ns*2, anc_nl, nd, value=-2, type=anc_type)
  ;seqpath = -1*ones(nl, ns*2);
  seqpath = make_array(anc_ns*2, anc_nl, value=-2, /integer)
  ;newseqpath = -1*ones(nl, ns*2);
  newseqpath = make_array(anc_ns*2, anc_nl, value=-2, /integer)
  
  ;seq1=zeros(nd, ns);
  seq1 = make_array(anc_ns, nd, type=anc_type)
  ;tmpind = true(1, ns);
  tmpind = bytarr(anc_ns)+1b
  ;for d=1:nd
  ;    seq1(d, :) = InputImage(2,:,d);
  ;    tmpind = tmpind & seq1(d,:)~=0;
  ;end
  for d=0,nd-1,1 do begin
    seq1[*,d] = wfmax_dual[*, 1, d]
    tmpind = tmpind and (seq1[*,d] ne 0)
  endfor
  ;seq1 = seq1(:, tmpind);
  seq1 = seq1[where(tmpind), *]
  ;
  ;seq2=zeros(nd, ns);
  seq2 = make_array(anc_ns, nd, type=anc_type)
  ;tmpind = true(1, ns);
  tmpind = bytarr(anc_ns)+1b
  ;for d=1:nd
  ;    seq2(d, :) = InputImage(3,:,d);
  ;    tmpind = tmpind & seq2(d,:)~=0;
  ;end
  for d=0,nd-1,1 do begin
    seq2[*,d] = wfmax_dual[*, 2, d]
    tmpind = tmpind and (seq2[*,d] ne 0)
  endfor
  ;seq2 = seq2(:, tmpind);
  seq2 = seq2[where(tmpind), *]
  ;
  ;[newseq1, newseq2, newseqlen, gfnewseq1, gfnewseq2, path1, path2] = lz_nwalign_nd(seq1, seq2, gapopen);
  tmpalign = numeric_nwalign_nd(seq1, seq2, gapopen)
  ;gfImage(2, 1:newseqlen, :) = shiftdim(newseq1', -1);
  ;gfImage(3, 1:newseqlen, :) = shiftdim(newseq2', -1);
  for d=0,nd-1,1 do begin
    gfImage[0:tmpalign.newseqlen-1, 1, d] = tmpalign.newseq1[*, d]
    gfImage[0:tmpalign.newseqlen-1, 2, d] = tmpalign.newseq2[*, d]
  endfor
  ;seqpath(2, 1:newseqlen) = path1;
  seqpath[0:tmpalign.newseqlen-1, 1] = tmpalign.path1
  ;seqpath(3, 1:newseqlen) = path2;
  seqpath[0:tmpalign.newseqlen-1, 2] = tmpalign.path2
  ;seq1 = gfnewseq2;
  seq1 = tmpalign.gfnewseq2
  ;for n=3:nl-1
  ;    seq2=zeros(nd, ns);
  ;    tmpind = true(1, ns);
  ;    for d=1:nd
  ;        seq2(d, :) = InputImage(n+1,:,d);
  ;        tmpind = tmpind & seq2(d,:)~=0;
  ;    end
  ;    seq2 = seq2(:, tmpind);
  ;
  ;    [newseq1, newseq2, newseqlen, gfnewseq1, gfnewseq2, path1, path2] = lz_nwalign_nd(seq1, seq2, gapopen);
  ;    newseqpath(n, 1:newseqlen) = path1;
  ;    seqpath(n+1, 1:newseqlen) = path2;
  ;    gfImage(n+1, 1:newseqlen, :) = shiftdim(newseq2', -1);
  ;    seq1 = gfnewseq2;
  ;    fprintf('%d\n', n);
  ;end
  n=2
  for n=2,anc_nl-2,1 do begin
    seq2 = make_array(anc_ns, nd, type=anc_type)
    tmpind = bytarr(anc_ns)+1b
    for d=0,nd-1,1 do begin
      seq2[*,d] = wfmax_dual[*, n+1, d]
      tmpind = tmpind and (seq2[*,d] ne 0)
    endfor
    seq2 = seq2[where(tmpind), *]
    tmpalign = numeric_nwalign_nd(seq1, seq2, gapopen)
    newseqpath[0:tmpalign.newseqlen-1, n] = tmpalign.path1
    seqpath[0:tmpalign.newseqlen-1, n+1] = tmpalign.path2
    for d=0,nd-1,1 do begin
      gfImage[0:tmpalign.newseqlen-1,n+1,d] = tmpalign.newseq2[*, d]
    endfor
    seq1 = tmpalign.gfnewseq2
    print, n
  endfor
  ;
  ;; update the path of the original sequence with the path of the gap-filled
  ;; sequence from last line to the first line
  ;AlignedImage = -1*ones(nl, ns*2, nd);
  AlignedImage = make_array(anc_ns*2, anc_nl, nd, value=-2, type=anc_type)
  ;AlignedImage(nl,:,:) = gfImage(nl,:,:);
  AlignedImage[*,anc_nl-1,*] = gfImage[*,anc_nl-1,*]
  ;
  ;p2 = newseqpath(nl-1,:); p2 = p2(p2~=-1);
  p2 = newseqpath[*,anc_nl-2]
  p2 = p2[where(p2 ne -2)]
  ;for n=nl-1:-1:4
  ;%     for d=1:nd
  ;%         tmpseq = gfImage(n,:,d); tmpseq = tmpseq(tmpseq~=-1);
  ;%         tmpseq2 = zeros(size(p2));
  ;%         tmpseq2(p2~=0) = tmpseq;
  ;%         AlignedImage(n, 1:length(p2),d) = shiftdim(tmpseq2', -1);
  ;%     end
  ;    tmpseq = gfImage(n,:,:);
  ;    tmpseq = tmpseq(1, gfImage(n,:,1)~=-1, :);
  ;    tmpseq2 = zeros(1, length(p2), nd);
  ;    tmpseq2(1, p2~=0, :) = tmpseq;
  ;    AlignedImage(n, 1:length(p2),:) = tmpseq2;
  ;
  ;    p1 = newseqpath(n-1, :); p1 = p1(p1~=-1);
  ;    tmpp1 = zeros(size(p2));
  ;    tmpp1(p2~=0) = p1;
  ;    p2 = tmpp1;
  ;end
  for n=anc_nl-2,3,-1 do begin
    for d=0,nd-1,1 do begin
      tmpseq = gfImage[*,n,d]
      tmpseq = tmpseq[where(tmpseq ne -2)]
      tmpseq2 = make_array(size(p2, /dimensions), type=anc_type)
      tmpseq2[where(p2 ne -1)] = tmpseq
      AlignedImage[0:size(p2, /n_elements)-1, n, d] = tmpseq2
    endfor
    p1 = newseqpath[*, n-1]
    p1 = p1[where(p1 ne -2)]
    tmpp1 = make_array(size(p2, /dimensions), /integer)
    tmpp1[where(p2 ne -1)] = p1
    p2 = tmpp1
  endfor
  
  ;for d=1:nd
  ;    tmpseq = gfImage(3,:,d); tmpseq = tmpseq(tmpseq~=-1);
  ;    tmpseq2 = zeros(size(p2));
  ;    tmpseq2(p2~=0) = tmpseq;
  ;    AlignedImage(3, 1:length(p2),d) = shiftdim(tmpseq2', -1);
  ;end
  ;for d=1:nd
  ;    tmpseq = gfImage(2,:,d); tmpseq = tmpseq(tmpseq~=-1);
  ;    tmpseq2 = zeros(size(p2));
  ;    tmpseq2(p2~=0) = tmpseq;
  ;    AlignedImage(2, 1:length(p2),d) = shiftdim(tmpseq2', -1);
  ;end
  for d=0,nd-1,1 do begin
    tmpseq = gfImage[*, 2, d]
    tmpseq = tmpseq[where(tmpseq ne -2)]
    tmpseq2 = make_array(size(p2, /dimensions), type=anc_type)
    tmpseq2[where(p2 ne -1)] = tmpseq
    AlignedImage[0:size(p2, /n_elements)-1, 2, d] = tmpseq2
  endfor
  for d=0,nd-1,1 do begin
    tmpseq = gfImage[*, 1, d]
    tmpseq = tmpseq[where(tmpseq ne -2)]
    tmpseq2 = make_array(size(p2, /dimensions), type=anc_type)
    tmpseq2[where(p2 ne -1)] = tmpseq
    AlignedImage[0:size(p2, /n_elements)-1, 1, d] = tmpseq2
  endfor
  ;
  ;AlignedImage = AlignedImage(:, 1:length(p2), :);
  AlignedImage = AlignedImage[0:size(p2, /n_elements)-1, *, *]
  ;tmpind = find(AlignedImage(2,:,1)~=0 & AlignedImage(2,:,1)~=-1);
  tmpind = where(AlignedImage[*,1,0] ne -1 and AlignedImage[*,1,0] ne -2)
  numind = size(tmpind, /n_elements)
  ;tmpcount = sum(InputImage(1, 1:ns, 1)~=0);
  tmpcount = total(wfmax_dual[0:anc_ns-1,0,0] ne 0)
  ;AlignedImage(1, tmpind(end-tmpcount+1:end), :) = InputImage(1,InputImage(1,1:ns,1)~=0,:);
  AlignedImage[tmpind[numind-tmpcount:numind-1], 0, *] = wfmax_dual[where(wfmax_dual[0:anc_ns-1,0,1] ne 0), 0, *]
  
  cleanup:
  
end