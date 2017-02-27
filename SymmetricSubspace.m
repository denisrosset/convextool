classdef SymmetricSubspace
    properties
        d;
        n;
        dim;
        cumProdFull;
    end
    methods
        function obj = SymmetricSubspace(d, n)
            obj.d = d;
            obj.n = n;
            obj.dim = SymmetricSubspace.computeDimension(d, n);
            p = cumprod(d*ones(1, n));
            obj.cumProdFull = [1 p(1:end-1)];
        end
        function sub = indToSubSym(obj, ind)
        % indToSubSym Convert from an index of the symmetric basis to subindices
        %
        % Let dim be the dimension of this SymmetricSubspace(d, n)
        %
        % INPUT  ind is a vector of m indices of the symmetric basis, values in 1...dim
        % OUTPUT sub is a (m x n) matrix, the r-th row is defined with respect to ind(r)
        %        sub(r, :) = [c1 c2 ... cn] with 1 <= c1 <= <= c2 <= ... <= cn <= d
            ind = ind(:);
            nRows = length(ind);
            if obj.n == 1
                sub = ind;
            else
                [sortedInd, index] = sort(ind);
                rStart = 1;
                nElementsBefore = 0;
                sortedSub = zeros(nRows, obj.n);
                for firstCol = 1:obj.d
                    % when the first column index is i, the elements of the remaining
                    % columns can be chosen from i to d, thus there are (d - i + 1)
                    % choices for these (n - 1) columns
                    sizeOfBlock = SymmetricSubspace.computeDimension(obj.d - firstCol + 1, obj.n - 1);
                    startIndex = nElementsBefore + 1;
                    endIndex = startIndex + sizeOfBlock - 1;
                    if rStart <= nRows && sortedInd(rStart) <= endIndex
                        % we have a block
                        rNextStart = rStart + 1;
                        while rNextStart <= nRows && sortedInd(rNextStart) <= endIndex
                            rNextStart = rNextStart + 1;
                        end
                        rEnd = rNextStart - 1;
                        remainingColsInd = sortedInd(rStart:rEnd) - nElementsBefore;
                        remainingColsSub = SymmetricSubspace(obj.d - firstCol + 1, obj.n - 1).indToSubSym(remainingColsInd);
                        sortedSub(rStart:rEnd, 1) = firstCol;
                        sortedSub(rStart:rEnd, 2:end) = remainingColsSub + firstCol - 1;
                        rStart = rNextStart;
                    end
                    nElementsBefore = nElementsBefore + sizeOfBlock;
                end
                sub = zeros(nRows, obj.n);
                sub(index, :) = sortedSub;
            end
        end
        function ind = subToIndSym(obj, sub)
        % subToIndSym Converts from subindices to the index in the symmetric basis
        %
        % Let dim be the dimension of this SymmetricSubspace(d, n)
        %
        % INPUT:  sub is a (m x n) matrix, whose r-th row is composed of subindices
        %         1 <= c1 ... cn <= d
        %
        % OUTPUT  ind is a vector of length m, whose r-th element correspond to the 
        %         index of sub(r, :) in the symmetric basis
            if obj.n == 1
                ind = sub;
            else
                sub = sort(sub')'; % put each row of subindices in increasing order
                                   % sort treats columns individually, so transpose
                [sortedSub, index] = sortrows(sub); % sorts the subindices for speed
                                                    % we have sortedSub = sub(index, :)
                rStart = 1;
                nRows = size(sortedSub, 1);
                sortedInd = zeros(nRows, 1);
                while rStart <= nRows
                    rNextStart = rStart + 1;
                    firstCol = sortedSub(rStart, 1);
                    while rNextStart <= nRows && sortedSub(rNextStart, 1) == firstCol
                        rNextStart = rNextStart + 1;
                    end
                    rEnd = rNextStart - 1;
                    nElementsBefore = 0;
                    for firstColBefore = 1:(firstCol-1)
                        % when the first column index is i, the elements of the remaining
                        % columns can be chosen from i to d, thus there are (d - i + 1)
                        % choices for these (n - 1) columns
                        sizeOfBlock = SymmetricSubspace.computeDimension(obj.d - firstColBefore + 1, obj.n - 1);
                        nElementsBefore = nElementsBefore + sizeOfBlock;
                    end
                    remainingCols = sortedSub(rStart:rEnd, 2:end) - firstCol + 1;
                    remainingColsInd = SymmetricSubspace(obj.d - firstCol + 1, obj.n - 1).subToIndSym(remainingCols);
                    sortedInd(rStart:rEnd) = nElementsBefore + remainingColsInd - 1;
                    rStart = rNextStart;
                end
                ind = zeros(nRows, 1);
                ind(index) = sortedInd + 1;
            end
        end
        function sub = indToSubFull(obj, ind)
            ind = ind(:) - 1;
            sub = zeros(length(ind), obj.n);
            for i = 1:obj.n
                sub(:, i) = mod(ind, obj.d);
                ind = ind - sub(:, i);
                sub(:, i) = sub(:, i) + 1;
                ind = ind ./ obj.d;
            end
        end
        function ind = subToIndFull(obj, sub);
            assert(size(sub, 2) == obj.n);
            ind = (obj.cumProdFull * (sub' - 1) + 1)';
        end
    end
    methods(Static)
        function dim = computeDimension(d, n)
            dim = nchoosek(n + d - 1, n);
        end
        function p = uniquePermutations(a)
        % Returns all the unique permutations of the vector a
        %
        % taken from https://www.mathworks.com/matlabcentral/newsreader/view_thread/164470
            [u, ~, J] = unique(a);
            p = u(up(J, length(a)));
            function p = up(J, n)
                ktab = histc(J,1:max(J));
                l = n;
                p = zeros(1, n);
                s = 1;
                for i=1:length(ktab)
                    k = ktab(i);
                    c = nchoosek(1:l, k);
                    m = size(c,1);
                    [t, ~] = find(~p.');
                    t = reshape(t, [], s);
                    c = t(c,:)';
                    s = s*m;
                    r = repmat((1:s)',[1 k]);
                    q = accumarray([r(:) c(:)], i, [s n]);
                    p = repmat(p, [m 1]) + q;
                    l = l - k;
                end
            end            
        end       
    end
end
