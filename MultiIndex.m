classdef MultiIndex
% MultiIndex provides efficient implementations of ind2sub/sub2ind
% by precomputing constants based on the dimensions, vectorizing
% operations and avoiding the use of cell arrays
    properties(SetAccess = immutable)
        n;
        dims;
        cumProd;
    end
    methods
        function obj = MultiIndex(dims)
            dims = dims(:)';
            obj.dims = dims;
            obj.n = length(dims);
            cp = cumprod(dims);
            obj.cumProd = [1 cp(1:end-1)];
        end
        function sub = indToSub(obj, ind)
        % given a vector of m indices, returns a matrix of size (m x n) of subindices
            ind = ind(:);
            nRows = length(ind);
            switch obj.n
              case 0
                sub = zeros(nRows, 0);
              case 1
                sub = ind;
              otherwise
                ind = ind - 1; % convert to 0-based indices
                sub = zeros(nRows, obj.n);
                for i = 1:obj.n
                    sub(:, i) = mod(ind, obj.dims(i));
                    ind = ind - sub(:, i);
                    ind = ind ./ obj.dims(i);
                end
                sub = sub + 1; % convert back to 1-based indices
            end
        end
        function ind = subToInd(obj, sub)
        % given a (m x n) matrix of subindices, returns the vector of the m corresponding indices
            nRows = size(sub, 1);
            switch obj.n
              case 0
                ind = ones(nRows, 1);
              case 1
                ind = sub;
              otherwise
                assert(size(sub, 2) == obj.n);
                ind = (obj.cumProd * (sub' - 1) + 1)';
            end
        end
    end
end
