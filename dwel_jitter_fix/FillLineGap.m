function gapline = FillLineGap(gapline, gappos)
% fill the gap in each sequence by linear interpolation from the closest
% neighbors. 
%gappos = find(path(:,1)==0);
seqlen = length(gapline);
if length(gappos)>=2
    tmp = gappos(2:end)-gappos(1:end-1);
    nseg = sum(tmp>1)+1;
    segpos = find(tmp>1)+1; segpos = reshape(segpos, 1, length(segpos));
    leftind = gappos([1, segpos]);
    rightind = gappos([segpos-1, length(gappos)]);
    if leftind(1)<=1
        gapline(:, leftind(1):rightind(1)) = repmat(gapline(:, rightind(1)+1), 1, rightind(1)-leftind(1)+1);
    elseif rightind(1)>=seqlen
        gapline(:, leftind(1):rightind(1)) = repmat(gapline(:, leftind(1)-1), 1, rightind(1)-leftind(1)+1);
    else
        gapline(:, leftind(1):rightind(1)) = repmat(0.5*(gapline(:, leftind(1)-1)+gapline(:, rightind(1)+1)), 1, rightind(1)-leftind(1)+1);
    end
    if nseg>=3
        for ns=2:nseg-1
            gapline(:, leftind(ns):rightind(ns)) = repmat(0.5*(gapline(:, leftind(ns)-1)+gapline(:, rightind(ns)+1)), 1, rightind(ns)-leftind(ns)+1);
        end
    end
    if nseg>=2
        if rightind(nseg)>=seqlen
            gapline(:, leftind(nseg):rightind(nseg)) = repmat(gapline(:, leftind(nseg)-1), 1, rightind(nseg)-leftind(nseg)+1);
        else
            gapline(:, leftind(nseg):rightind(nseg)) = repmat(0.5*(gapline(:, leftind(nseg)-1)+gapline(:, rightind(nseg)+1)), 1, rightind(nseg)-leftind(nseg)+1);
        end
    end
else
    if length(gappos)>=1
        if gappos(1)==1
            gapline(:, gappos(1)) = gapline(:, gappos(1)+1);
        elseif gappos(1)==seqlen
            gapline(:, gappos(1)) = gapline(:, gappos(1)-1);
        else
            gapline(:, gappos(1)) = 0.5*(gapline(:, gappos(1)-1)+gapline(:, gappos(1)+1));
        end
    end
end
end