function simplegap_nd, seq1, m, seq2, n, gap
  ; Variant of Standard Needleman-Wunsch algorithm by Zhan Li for scan line
  ; alignment.
  
  ; set up storage for dynamic programming matrix
  ;F = zeros(n+1,m+1)
  F = make_array(m+1, n+1, /float)
  ;F(2:end,1) = gap * (1:n)';
  F[0, 1:n] = gap * transpose(indgen(n)+1)
  ;F(1,2:end) = gap * (1:m);
  F[1:m, 0] = gap * (indgen(m)+1)
  ;
  ;; and for the back tracing matrix
  ;pointer= repmat(uint8(4),n+1,m+1);
  pointer = make_array(m+1, n+1, /byte, value=4)
  ;pointer(:,1) = 2;  % up
  pointer[0,*] = 2
  ;pointer(1,1) = 1;  % diagonal
  pointer[0,0] = 1
  ;
  ;; initialize buffers to the first column
  ;ptr = pointer(:,2); % ptr(1) is always 4
  ptr = pointer[1, *] ; ptr[0] is always 4
  ;currentFColumn = F(:,1);
  currentFColumn = F[0,*]
  ;
  ;; main loop runs through the matrix looking for maximal scores
  ;for outer = 2:m+1
  for outer = 1, m, 1 do begin
    ;    ; score current column
    ;    %scoredMatchColumn = ScoringMatrix(intseq2,intseq1(outer-1));
    ;    scoredMatchColumn = -1*sum((intseq2 - repmat(intseq1(:, outer-1), 1, n)).^2, 1);
    scoredMatchColumn = -1*total((seq2 - transpose(cmreplicate(transpose(seq1[outer-1, *]), [n, 1])))^2, 2)
    ;    ; grab the data from the matrices and initialize some values
    ;    lastFColumn    = currentFColumn;
    lastFColumn = currentFColumn
    ;    currentFColumn = F(:,outer);
    currentFColumn = F[outer, *]
    ;    best = currentFColumn(1);
    best = currentFColumn[0]
    
    ;    for inner = 2:n+1
    for inner = 1,n,1 do begin
      ;        ; score the three options
      ;        up       = best + gap;
      up = best + gap
      ;        left     = lastFColumn(inner) + gap;
      left = lastFColumn[inner]+gap
      ;        diagonal = lastFColumn(inner-1) + scoredMatchColumn(inner-1);
      diagonal = lastFColumn[inner-1] + scoredMatchColumn[inner-1]
      
      ;        ; max could be used here but it is quicker to use if statements
      ;        if up > left
      ;            best = up;
      ;            pos = 2;
      ;        else
      ;            best = left;
      ;            pos = 4;
      ;        end
      if up gt left then begin
        best = up
        pos = 2
      endif else begin
        best = left
        pos = 4
      endelse
      
      ;        if diagonal >= best
      ;            best = diagonal;
      ;            ptr(inner) = 1;
      ;        else
      ;            ptr(inner) = pos;
      ;        end
      if diagonal ge best then begin
        best = diagonal
        ptr[inner] = 1
      endif else begin
        ptr[inner] = pos
      endelse
      ;        currentFColumn(inner) = best;
      currentFColumn[inner] = best
      
      ;    end % inner
    endfor ; inner
    ;    ; put back updated columns
    ;    F(:,outer)   = currentFColumn;
    F[outer, *] = currentFColumn
    ;    ; save columns of pointers
    pointer[outer, *] = ptr
    ;    pointer(:,outer)  = ptr;
    ;end % outer
  endfor ;outer
  
  ;; get the alignment path from the pointer
  ;path = zeros(n+m, 2);
  path = make_array(2, n+m, /integer, value=-1)
  ;mcol = m+1;
  mcol = m
  ;nrow = n+1;
  nrow = n
  
  ;for p=1:n+m
  for p = 0,n+m-1,1 do begin
    ;    if (mcol<1 || nrow<1)
    ;        break;
    ;    end
    if (mcol lt 1) || (nrow lt 1) then begin
      break
    endif
    ;    if (pointer(nrow, mcol)==1) % diagonal, aligned
    ;        path(p, :) = [mcol-1, nrow-1];
    ;        mcol = mcol - 1;
    ;        nrow = nrow - 1;
    ;    elseif (pointer(nrow, mcol)==2) % up, gap in seq1 (m)
    ;        path(p, :) = [0, nrow-1];
    ;        nrow = nrow - 1;
    ;    else % left, gap in seq 2 (n)
    ;        path(p, :) = [mcol-1, 0];
    ;        mcol = mcol - 1;
    ;    end
    if pointer[mcol, nrow] eq 1 then begin ; diagonal, aligned
      path[*,p] = [mcol-1, nrow-1]
      mcol = mcol - 1
      nrow = nrow - 1
    endif else begin
      if pointer[mcol, nrow] eq 2 then begin ; up, gap in seq1 (m)
        path[*,p] = [-1, nrow-1]
        nrow = nrow - 1
      endif else begin; left, gap in seq2 (n)
        path[*,p] = [mcol-1, -1]
        mcol = mcol - 1
      endelse
    endelse
    ;end % p
  endfor ; p
  return, {F:F, pointer:pointer, path:path}
end

function FillLineGap, gapline, gappos
  ; fill the gap in each sequence by linear interpolation from the closest
  ; neighbors.
  ; gapline and gappos are both row vectors
  
  ;seqlen = length(gapline);
  seqlen = size(gapline, /n_elements)
  gapnum = size(gappos, /n_elements)
  ;if ~isempty(gappos)
  if gapnum ge 2 then begin
    ;    tmp = gappos(2:end)-gappos(1:end-1);
    tmp = gappos[1:gapnum-1] - gappos[0:gapnum-2]
    ;    nseg = sum(tmp>1)+1;
    ;    segpos = find(tmp>1)+1; segpos = reshape(segpos, 1, length(segpos));
    segpos = where(tmp gt 1, nseg) + 1
    nseg = nseg + 1
    ;    leftind = gappos([1, segpos]);
    leftind = gappos[[0, segpos]]
    ;    rightind = gappos([segpos-1, length(gappos)]);
    rightind = gappos[[segpos-1, gapnum]]
    ;    if leftind(1)<=1
    ;        gapline(:, leftind(1):rightind(1)) = repmat(gapline(:, rightind(1)+1), 1, rightind(1)-leftind(1)+1);
    ;    elseif rightind(1)>=seqlen
    ;        gapline(:, leftind(1):rightind(1)) = repmat(gapline(:, leftind(1)-1), 1, rightind(1)-leftind(1)+1);
    ;    else
    ;        gapline(:, leftind(1):rightind(1)) = repmat(0.5*(gapline(:, leftind(1)-1)+gapline(:, rightind(1)+1)), 1, rightind(1)-leftind(1)+1);
    ;    end
    if leftind[0] le 0 then begin
      gapline[leftind[0]:rightind[0]] = gapline[rightind[0]+1]
    endif else begin
      if rightind[0] ge seqlen-1 then begin
        gapline[leftind[0]:rightind[0]] = gapline[leftind[0]-1]
      endif else begin
        gapline[leftind[0]:rightind[0]] = 0.5*(gapline[leftind[0]-1]+gapline[rightind[0]+1])
      endelse
    endelse
    
    ;    if nseg>=3
    ;        for ns=2:nseg-1
    ;            gapline(:, leftind(ns):rightind(ns)) = repmat(0.5*(gapline(:, leftind(ns)-1)+gapline(:, rightind(ns)+1)), 1, rightind(ns)-leftind(ns)+1);
    ;        end
    ;    end
    ;    if nseg>=2
    ;        if rightind(nseg)>=seqlen
    ;            gapline(:, leftind(nseg):rightind(nseg)) = repmat(gapline(:, leftind(nseg)-1), 1, rightind(nseg)-leftind(nseg)+1);
    ;        else
    ;            gapline(:, leftind(nseg):rightind(nseg)) = repmat(0.5*(gapline(:, leftind(nseg)-1)+gapline(:, rightind(nseg)+1)), 1, rightind(nseg)-leftind(nseg)+1);
    ;        end
    ;    end
    if nseg ge 3 then begin
      for ns=1,nseg-2,1 do begin
        gapline[leftind[ns]:rightind[ns]] = 0.5*(gapline[leftind[ns]-1]+gapline[rightind[ns]+1])
      endfor
    endif
    if nseg ge 2 then begin
      if rightind[nseg-1] ge seqlen-1 then begin
        gapline[leftind[nseg-1]:rightind[nseg-1]] = gapline[leftind[nseg-1]-1]
      endif else begin
        gapline[leftind[nseg-1]:rightind[nseg-1]] = 0.5*(gapline[leftind[nseg-1]-1]+gapline[rightind[nseg-1]+1])
      endelse
    endif
  endif else begin
    if gapnum ge 1 then begin
      if gappos[0] eq 0 then begin
        gapline[gappos[0]] = gapline[gappos[0]+1]
      endif else begin
        if gappos[0] eq seqlen-1 then begin
          gapline[gappos[0]] = gapline[gappos[0]-1]
        endif else begin
          gapline[gappos[0]] = 0.5*(gapline[gappos[0]-1] + gapline[gappos[0]+1])
        endelse
      endelse
    endif
  endelse
  return, gapline
end

function numeric_nwalign_nd, seq1, seq2, gapopen
  ;NWALIGN performs Needleman-Wunsch global alignment of two sequences in
  ;numeric data type for DWEL's scan line alignment
  ;   References:
  ;   R. Durbin, S. Eddy, A. Krogh, and G. Mitchison. Biological Sequence
  ;   Analysis. Cambridge UP, 1998.
  ;   Needleman, S. B., Wunsch, C. D., J. Mol. Biol. (1970) 48:443-453
  ;INPUT:
  ; seq1: d*m, d is the dimension of the sequence, m is the length of the
  ; sequence
  ; seq2: d*n, d is the dimension of the sequence, n is the length of the
  ; sequence.
  
  tmp = size(seq1, /dimensions)
  m = tmp[0]
  nd = tmp[1]
  tmp = size(seq2, /dimensions)
  n = tmp[0]
  if (m eq 0) || (n eq 0) then begin
    print, 'numeric_nwalign_nd: Length of input sequences must be greater than 0.'
    return, -1
  endif
  if nd ne tmp[1] then begin
    print, 'numeric_nwalign_nd: Dimensions of two input sequences must be the same.'
    return, -1
  endif
  
  ;[F, pointer, path]=simplegap2(intseq1,m,intseq2,n,gapopen);
  nwalign = simplegap_nd(seq1, m, seq2, n, gapopen)
  path = nwalign.path
  
  ;path = path(sum(path,2)>0,:);
  path = path[*, where(total(path, 1) gt -2)]
  ;path = flipud(path);
  path = reverse(path, 2)
  
  ;newseqlen = size(path, 1);
  tmp = size(path, /dimensions)
  newseqlen = tmp[1]
  
  ;path1 = reshape(path(:,1), 1, newseqlen);
  ;path2 = reshape(path(:,2), 1, newseqlen);
  path1 = transpose(path[0,*])
  path2 = transpose(path[1,*])
  
  ;newseq1 = zeros(d, newseqlen);
  ;newseq2 = newseq1;
  seqtype = size(seq1, /type)
  newseq1 = make_array(newseqlen, nd, type=seqtype)
  newseq2 = newseq1
  ;newseq1(:, path(:,1)~=0) = seq1(:, path(path(:,1)~=0, 1));
  ;newseq2(:, path(:,2)~=0) = seq2(:, path(path(:,2)~=0, 2));
  newseq1[where(path1 ne -1), *] = seq1[path1[where(path1 ne -1)], *]
  newseq2[where(path2 ne -1), *] = seq2[path2[where(path2 ne -1)], *]
  
  ;gfnewseq1 = newseq1;
  ;gfnewseq2 = newseq2;
  gfnewseq1 = newseq1
  gfnewseq2 = newseq2
  
  ;; fill the gap in each sequence by linear interpolation from the closest
  ;; neighbors.
  ;gappos = find(path(:,1)==0);
  gappos = where(path1 eq -1)
  for d=0,nd-1,1 do begin
    gfnewseq1[*,d] = FillLineGap(gfnewseq1[*, d], gappos)
  endfor
  
  ;if ~isempty(gappos)
  ;    tmp = gappos(2:end)-gappos(1:end-1);
  ;    nseg = sum(tmp>1)+1;
  ;    segpos = find(tmp>1)+1; segpos = reshape(segpos, 1, length(segpos));
  ;    leftind = gappos([1, segpos]);
  ;    rightind = gappos([segpos-1, length(gappos)]);
  ;    if leftind(1)<=1
  ;        gfnewseq1(:, leftind(1):rightind(1)) = repmat(gfnewseq1(:, rightind(1)+1), 1, rightind(1)-leftind(1)+1);
  ;    elseif rightind(1)>=newseqlen
  ;        gfnewseq1(:, leftind(1):rightind(1)) = repmat(gfnewseq1(:, leftind(1)-1), 1, rightind(1)-leftind(1)+1);
  ;    else
  ;        gfnewseq1(:, leftind(1):rightind(1)) = repmat(0.5*(gfnewseq1(:, leftind(1)-1)+gfnewseq1(:, rightind(1)+1)), 1, rightind(1)-leftind(1)+1);
  ;    end
  ;    if nseg>=3
  ;        for ns=2:nseg-1
  ;            gfnewseq1(:, leftind(ns):rightind(ns)) = repmat(0.5*(gfnewseq1(:, leftind(ns)-1)+gfnewseq1(:, rightind(ns)+1)), 1, rightind(ns)-leftind(ns)+1);
  ;        end
  ;    end
  ;    if nseg>=2
  ;        if rightind(nseg)>=newseqlen
  ;            gfnewseq1(:, leftind(nseg):rightind(nseg)) = repmat(gfnewseq1(:, leftind(nseg)-1), 1, rightind(nseg)-leftind(nseg)+1);
  ;        else
  ;            gfnewseq1(:, leftind(nseg):rightind(nseg)) = repmat(0.5*(gfnewseq1(:, leftind(nseg)-1)+gfnewseq1(:, rightind(nseg)+1)), 1, rightind(nseg)-leftind(nseg)+1);
  ;        end
  ;    end
  ;end
  ;
  ;gappos = find(path(:,2)==0);
  gappos = where(path2 eq 0)
  for d=0,nd-1,1 do begin
    gfnewseq2[*,d] = FillLineGap(gfnewseq2[*, d], gappos)
  endfor
  ;if ~isempty(gappos)
  ;    tmp = gappos(2:end)-gappos(1:end-1);
  ;    nseg = sum(tmp>1)+1;
  ;    segpos = find(tmp>1)+1; segpos = reshape(segpos, 1, length(segpos));
  ;    leftind = gappos([1, segpos]);
  ;    rightind = gappos([segpos-1, length(gappos)]);
  ;    if leftind(1)<=1
  ;        gfnewseq2(:, leftind(1):rightind(1)) = repmat(gfnewseq2(:, rightind(1)+1), 1, rightind(1)-leftind(1)+1);
  ;    elseif rightind(1)>=newseqlen
  ;        gfnewseq2(:, leftind(1):rightind(1)) = repmat(gfnewseq2(:, leftind(1)-1), 1, rightind(1)-leftind(1)+1);
  ;    else
  ;        gfnewseq2(:, leftind(1):rightind(1)) = repmat(0.5*(gfnewseq2(:, leftind(1)-1)+gfnewseq2(:, rightind(1)+1)), 1, rightind(1)-leftind(1)+1);
  ;    end
  ;    if nseg>=3
  ;        for ns=2:nseg-1
  ;            gfnewseq2(:, leftind(ns):rightind(ns)) = repmat(0.5*(gfnewseq2(:, leftind(ns)-1)+gfnewseq2(:, rightind(ns)+1)), 1, rightind(ns)-leftind(ns)+1);
  ;        end
  ;    end
  ;    if nseg>=2
  ;        if rightind(nseg)>=newseqlen
  ;            gfnewseq2(:, leftind(nseg):rightind(nseg)) = repmat(gfnewseq2(:, leftind(nseg)-1), 1, rightind(nseg)-leftind(nseg)+1);
  ;        else
  ;            gfnewseq2(:, leftind(nseg):rightind(nseg)) = repmat(0.5*(gfnewseq2(:, leftind(nseg)-1)+gfnewseq2(:, rightind(nseg)+1)), 1, rightind(nseg)-leftind(nseg)+1);
  ;        end
  ;    end
  ;end
  return, {newseq1:newseq1, newseq2:newseq2, newseqlen:newseqlen, gfnewseq1:gfnewseq1, gfnewseq2:gfnewseq2, path1:path1, path2:path2}
end