function [newseq1, newseq2, newseqlen, gfnewseq1, gfnewseq2, path1, path2] = lz_nwalign_nd(seq1,seq2,gapopen)
%NWALIGN performs Needleman-Wunsch global alignment of two sequences in
%numeric data type for DWEL's scan line alignment
%   References:
%   R. Durbin, S. Eddy, A. Krogh, and G. Mitchison. Biological Sequence
%   Analysis. Cambridge UP, 1998.
%   Needleman, S. B., Wunsch, C. D., J. Mol. Biol. (1970) 48:443-453
%INPUT:
% seq1: d*m, d is the dimension of the sequence, m is the length of the
% sequence
% seq2: d*n, d is the dimension of the sequence, n is the length of the
% sequence. 

m = size(seq1, 2);
n = size(seq2, 2);
if ~n||~m
    error('lz_nwalign2:InvalidLengthSequences','Length of input sequences must be greater than 0');
end
if size(seq1, 1)~=size(seq2, 1)
    error('lz_nwalign2:InvalidLengthSequences','Dimensions of two input sequences must be the same');
end

d = size(seq1, 1);

intseq1 = seq1;
intseq2 = seq2;

[F, pointer, path]=simplegap2(intseq1,m,intseq2,n,gapopen);

path = path(sum(path,2)>0,:);
path = flipud(path);

newseqlen = size(path, 1);

path1 = reshape(path(:,1), 1, newseqlen);
path2 = reshape(path(:,2), 1, newseqlen);

newseq1 = zeros(d, newseqlen);
newseq2 = newseq1;
newseq1(:, path(:,1)~=0) = seq1(:, path(path(:,1)~=0, 1));
newseq2(:, path(:,2)~=0) = seq2(:, path(path(:,2)~=0, 2));

gfnewseq1 = newseq1;
gfnewseq2 = newseq2;

% fill the gap in each sequence by linear interpolation from the closest
% neighbors. 
gappos = find(path(:,1)==0);
% if ~isempty(gappos)
%     tmp = gappos(2:end)-gappos(1:end-1);
%     nseg = sum(tmp>1)+1;
%     segpos = find(tmp>1)+1; segpos = reshape(segpos, 1, length(segpos));
%     leftind = gappos([1, segpos]);
%     rightind = gappos([segpos-1, length(gappos)]);
%     if leftind(1)<=1
%         gfnewseq1(:, leftind(1):rightind(1)) = repmat(gfnewseq1(:, rightind(1)+1), 1, rightind(1)-leftind(1)+1);
%     elseif rightind(1)>=newseqlen
%         gfnewseq1(:, leftind(1):rightind(1)) = repmat(gfnewseq1(:, leftind(1)-1), 1, rightind(1)-leftind(1)+1);
%     else
%         gfnewseq1(:, leftind(1):rightind(1)) = repmat(0.5*(gfnewseq1(:, leftind(1)-1)+gfnewseq1(:, rightind(1)+1)), 1, rightind(1)-leftind(1)+1);
%     end
%     if nseg>=3
%         for ns=2:nseg-1
%             gfnewseq1(:, leftind(ns):rightind(ns)) = repmat(0.5*(gfnewseq1(:, leftind(ns)-1)+gfnewseq1(:, rightind(ns)+1)), 1, rightind(ns)-leftind(ns)+1);
%         end
%     end
%     if nseg>=2
%         if rightind(nseg)>=newseqlen
%             gfnewseq1(:, leftind(nseg):rightind(nseg)) = repmat(gfnewseq1(:, leftind(nseg)-1), 1, rightind(nseg)-leftind(nseg)+1);
%         else
%             gfnewseq1(:, leftind(nseg):rightind(nseg)) = repmat(0.5*(gfnewseq1(:, leftind(nseg)-1)+gfnewseq1(:, rightind(nseg)+1)), 1, rightind(nseg)-leftind(nseg)+1);
%         end
%     end
% end
gfnewseq1 = FillLineGap(gfnewseq1, gappos);

gappos = find(path(:,2)==0);
% if ~isempty(gappos)
%     tmp = gappos(2:end)-gappos(1:end-1);
%     nseg = sum(tmp>1)+1;
%     segpos = find(tmp>1)+1; segpos = reshape(segpos, 1, length(segpos));
%     leftind = gappos([1, segpos]);
%     rightind = gappos([segpos-1, length(gappos)]);
%     if leftind(1)<=1
%         gfnewseq2(:, leftind(1):rightind(1)) = repmat(gfnewseq2(:, rightind(1)+1), 1, rightind(1)-leftind(1)+1);
%     elseif rightind(1)>=newseqlen
%         gfnewseq2(:, leftind(1):rightind(1)) = repmat(gfnewseq2(:, leftind(1)-1), 1, rightind(1)-leftind(1)+1);
%     else
%         gfnewseq2(:, leftind(1):rightind(1)) = repmat(0.5*(gfnewseq2(:, leftind(1)-1)+gfnewseq2(:, rightind(1)+1)), 1, rightind(1)-leftind(1)+1);
%     end
%     if nseg>=3
%         for ns=2:nseg-1
%             gfnewseq2(:, leftind(ns):rightind(ns)) = repmat(0.5*(gfnewseq2(:, leftind(ns)-1)+gfnewseq2(:, rightind(ns)+1)), 1, rightind(ns)-leftind(ns)+1);
%         end
%     end
%     if nseg>=2
%         if rightind(nseg)>=newseqlen
%             gfnewseq2(:, leftind(nseg):rightind(nseg)) = repmat(gfnewseq2(:, leftind(nseg)-1), 1, rightind(nseg)-leftind(nseg)+1);
%         else
%             gfnewseq2(:, leftind(nseg):rightind(nseg)) = repmat(0.5*(gfnewseq2(:, leftind(nseg)-1)+gfnewseq2(:, rightind(nseg)+1)), 1, rightind(nseg)-leftind(nseg)+1);
%         end
%     end
% end
gfnewseq2 = FillLineGap(gfnewseq2, gappos);


%=== SIMPLEGAP is now a mex function ===%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F, pointer, path] = simplegap2(intseq1,m,intseq2,n,gap)
% Variant of Standard Needleman-Wunsch algorithm by Zhan Li for scan line
% alignment. 

% set up storage for dynamic programming matrix
F = zeros(n+1,m+1);
F(2:end,1) = gap * (1:n)';
F(1,2:end) = gap * (1:m);
% F(2:end,1) = ones(n,1)*-Inf;
% F(1,2:end) = ones(1,m)*-Inf;

% and for the back tracing matrix
pointer= repmat(uint8(4),n+1,m+1);
pointer(:,1) = 2;  % up
pointer(1,1) = 1;  % diagonal

% initialize buffers to the first column
ptr = pointer(:,2); % ptr(1) is always 4
currentFColumn = F(:,1);

% main loop runs through the matrix looking for maximal scores
for outer = 2:m+1

    % score current column
    %scoredMatchColumn = ScoringMatrix(intseq2,intseq1(outer-1));
    scoredMatchColumn = -1*sum((intseq2 - repmat(intseq1(:, outer-1), 1, n)).^2, 1);
    % grab the data from the matrices and initialize some values
    lastFColumn    = currentFColumn;
    currentFColumn = F(:,outer);
    best = currentFColumn(1);

    for inner = 2:n+1
        % score the three options
        up       = best + gap;
        left     = lastFColumn(inner) + gap;
        diagonal = lastFColumn(inner-1) + scoredMatchColumn(inner-1);

        % max could be used here but it is quicker to use if statements
        if up > left
            best = up;
            pos = 2;
        else
            best = left;
            pos = 4;
        end

        if diagonal >= best
            best = diagonal;
            ptr(inner) = 1;
        else
            ptr(inner) = pos;
        end
        currentFColumn(inner) = best;

    end % inner
    % put back updated columns
    F(:,outer)   = currentFColumn;
    % save columns of pointers
    pointer(:,outer)  = ptr;
end % outer

% get the alignment path from the pointer
path = zeros(n+m, 2);
mcol = m+1;
nrow = n+1;

for p=1:n+m
    if (mcol<1 || nrow<1)
        break;
    end
    if (pointer(nrow, mcol)==1) % diagonal, aligned
        path(p, :) = [mcol-1, nrow-1];
        mcol = mcol - 1;
        nrow = nrow - 1;
    elseif (pointer(nrow, mcol)==2) % up, gap in seq1 (m)
        path(p, :) = [0, nrow-1];
        nrow = nrow - 1;
    else % left, gap in seq 2 (n)
        path(p, :) = [mcol-1, 0];
        mcol = mcol - 1;
    end
end % p
