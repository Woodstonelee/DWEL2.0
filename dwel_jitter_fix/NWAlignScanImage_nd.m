function [AlignedImage, AlignedMask] = NWAlignScanImage_nd(InputImage, gapopen, startline)
nl = size(InputImage, 1);
ns = size(InputImage, 2);
nd = size(InputImage, 3);

% use new line2 as the next line1 and do nothing to combine two versions of
% alignments. 
gfImage = -1*ones(nl, ns*2, nd);
seqpath = -1*ones(nl, ns*2);
newseqpath = -1*ones(nl, ns*2);

seq1=zeros(nd, ns);
tmpind = true(1, ns);
for d=1:nd
    seq1(d, :) = InputImage(startline,:,d);
    tmpind = tmpind & seq1(d,:)~=0;
end
seq1 = seq1(:, tmpind);

seq2=zeros(nd, ns);
tmpind = true(1, ns);
for d=1:nd
    seq2(d, :) = InputImage(startline+1,:,d);
    tmpind = tmpind & seq2(d,:)~=0;
end
seq2 = seq2(:, tmpind);

[newseq1, newseq2, newseqlen, gfnewseq1, gfnewseq2, path1, path2] = lz_nwalign_nd(seq1, seq2, gapopen);
gfImage(startline, 1:newseqlen, :) = shiftdim(newseq1', -1);
gfImage(startline+1, 1:newseqlen, :) = shiftdim(newseq2', -1);
seqpath(startline, 1:newseqlen) = path1;
seqpath(startline+1, 1:newseqlen) = path2;
seq1 = gfnewseq2;
strCR = [];
for n=startline+1:nl-1
    seq2=zeros(nd, ns);
    tmpind = true(1, ns);
    for d=1:nd
        seq2(d, :) = InputImage(n+1,:,d);
        tmpind = tmpind & seq2(d,:)~=0;
    end
    seq2 = seq2(:, tmpind);
    
    [newseq1, newseq2, newseqlen, gfnewseq1, gfnewseq2, path1, path2] = lz_nwalign_nd(seq1, seq2, gapopen);
    newseqpath(n, 1:newseqlen) = path1;
    seqpath(n+1, 1:newseqlen) = path2;
    gfImage(n+1, 1:newseqlen, :) = shiftdim(newseq2', -1);
    seq1 = gfnewseq2;
    strout = num2str(n);
    fprintf([strCR, strout]);
    strCR = repmat('\b', 1, length(strout));
end
fprintf('\n');

% update the path of the original sequence with the path of the gap-filled
% sequence from last line to the first line
AlignedImage = -1*ones(nl, ns*2, nd);
AlignedImage(nl,:,:) = gfImage(nl,:,:);

p2 = newseqpath(nl-1,:); p2 = p2(p2~=-1);
for n=nl-1:-1:startline+2
    tmpseq = gfImage(n,:,:); 
    tmpseq = tmpseq(1, gfImage(n,:,1)~=-1, :);
    tmpseq2 = zeros(1, length(p2), nd);
    tmpseq2(1, p2~=0, :) = tmpseq;
    AlignedImage(n, 1:length(p2),:) = tmpseq2;
    
    p1 = newseqpath(n-1, :); p1 = p1(p1~=-1);
    tmpp1 = zeros(size(p2));
    tmpp1(p2~=0) = p1;
    p2 = tmpp1;
end

for d=1:nd
    tmpseq = gfImage(startline+1,:,d); tmpseq = tmpseq(tmpseq~=-1);
    tmpseq2 = zeros(size(p2)); 
    tmpseq2(p2~=0) = tmpseq;
    AlignedImage(startline+1, 1:length(p2),d) = shiftdim(tmpseq2', -1);
end
for d=1:nd
    tmpseq = gfImage(startline,:,d); tmpseq = tmpseq(tmpseq~=-1);
    tmpseq2 = zeros(size(p2)); 
    tmpseq2(p2~=0) = tmpseq;
    AlignedImage(startline, 1:length(p2),d) = shiftdim(tmpseq2', -1);
end

AlignedImage = AlignedImage(:, 1:length(p2), :);
align_ns = length(p2);
if startline>1
    for il = 1:startline-1
        tmpind = InputImage(il, :, 1)~=0;
        if sum(tmpind(:))<=align_ns
            AlignedImage(il, end-sum(tmpind(:))+1:end, :) = InputImage(il, tmpind, :);
        else
            AlignedImage(il, :, :) = InputImage(il, tmpind(end-align_ns+1:end), :);
        end
    end 
end

AlignedMask = zeros(size(AlignedImage(:,:,1)));
AlignedMask(AlignedImage(:,:,1)~=0 & AlignedImage(:,:,1)~=-1) = 1;
end