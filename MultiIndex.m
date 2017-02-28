classdef MultiIndex
% MultiIndex provides efficient implementations of ind2sub/sub2ind
% vectorizing operations and avoiding the use of cell arrays
%
% Setting optimize = 1 will switch the algorithm to 
% "fast integer by precomputation of a reciprocal"; we lift the algorithm
% of libdivide.com
%
% For fastest results, call the functions with integer types such as 'uint32'
% or 'uint64'.
    properties(SetAccess = immutable)
        n;
        dims;
    end
    properties(Access = private)
        cumProd;
        % constants for fast integer division
        optimize;
        pre_shift;
        multiplier;
        increment;
        full_shift;
    end
    methods
        function obj = MultiIndex(dims, optimize)
            if nargin > 1
                obj.optimize = optimize;
            else
                obj.optimize = 0;
            end
            dims = dims(:)';
            obj.dims = dims;
            obj.n = length(dims);
            cp = cumprod(uint32(dims));
            obj.cumProd = [1 cp(1:end-1)];
            if obj.optimize && prod(dims) < 2^32 % everything fits in 32 bits
                obj.pre_shift = zeros(1, obj.n, 'int32');
                obj.multiplier = zeros(1, obj.n, 'uint64');
                obj.increment = zeros(1, obj.n, 'int32');
                obj.full_shift = zeros(1, obj.n, 'int32');
                for i = 1:obj.n
                    ld = LibDivide32(obj.dims(i));
                    obj.pre_shift(i) = ld.pre_shift;
                    obj.multiplier(i) = ld.multiplier;
                    obj.increment(i) = ld.increment;
                    obj.full_shift(i) = ld.full_shift;
                end
            end
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
                if ~obj.optimize
                    ind = uint32(ind - 1); % convert to 0-based indices
                    sub = zeros(nRows, obj.n, 'uint32');
                    for i = 1:obj.n
                        % TODO: replace modulo by libdivide-like method
                        sub(:, i) = mod(ind, uint32(obj.dims(i)));
                        ind = ind - sub(:, i);
                        ind = ind ./ obj.dims(i);
                    end
                    sub = cast(sub + 1, class(ind)); % convert back to 1-based indices
                else
                    sub = MultiIndex.indToSubOptimized(uint32(ind), ...
                                                       uint32(obj.dims), ...
                                                       obj.pre_shift, ...
                                                       obj.multiplier, ...
                                                       obj.increment, ...
                                                       obj.full_shift);
                end
            end
            sub = cast(sub, class(ind));
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
                nRows = size(sub, 1);
                assert(size(sub, 2) == obj.n);
                if isa(sub, 'uint32')
                    ind = zeros(nRows, 1, 'uint32');
                    ind = sub(:, 1);
                    for i = 2:n
                        ind = ind + obj.cumProd(i) * (sub(:, i) - 1);
                    end
                else
                    ind = (double(obj.cumProd) * (sub' - 1) + 1)';
                end
            end
        end
    end
    methods(Static)
        function sub = indToSubOptimized(ind, dims, pre_shift, multiplier, increment, full_shift)
            ind = ind - 1; % convert to 0-based indices
            n = length(dims);
            sub = zeros(length(ind), n, 'uint32');
            for i = 1:length(pre_shift)
                % we use the second "emitted code" suggestion from libdivide.c
                rem = ind;
                if pre_shift(i) > 0
                    rem = bitshift(rem, -pre_shift(i));
                end
                if increment(i)
                    rem1 = uint32(bitshift(multiplier(i) * uint64(rem) + multiplier(i), -full_shift(i)));
                else
                    rem1 = uint32(bitshift(multiplier(i) * uint64(rem), -full_shift(i)));
                end
                sub(:, i) = (ind - dims(i)*rem1) + 1;
                ind = rem;
            end
        end
    end
end
