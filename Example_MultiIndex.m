dims = [3 3];
ind = (1:9)'; % MultiIndex expects indices as column vectors

% standard Matlab function
[I J] = ind2sub(dims, ind);

% but when the number of dimensions is not fixed, it becomes crazy
out = cell(1, 2);
[out{:}] = ind2sub(dims, ind);
I = out{1};
J = out{2};

% faster indToSub method that returns an index array 
IJ = MultiIndex([3 3]).indToSub(ind);

% both return the same results in the end
assert(isequal(I, IJ(:, 1)));
assert(isequal(J, IJ(:, 2)));

% the reverse transformation

ind1 = sub2ind(dims, I, J);

% or, when the number of dimensions is not fixed, cell array again!
in = {I J};
ind2 = sub2ind(dims, in{:});

% while MultiIndex is more efficient/practical here
ind3 = MultiIndex(dims).subToInd(IJ);

% results are identical
assert(isequal(ind, ind1));
assert(isequal(ind, ind2));
assert(isequal(ind, ind3));

% if you do many small transformation with the same dimension, 
% think of creating MultiIndex only once

mi = MultiIndex(dims);
for i = 1:3
    for j = 1:3
        % ...
        mi.subToInd([i j]);
        % instead of
        MultiIndex(dims).subToInd([i j]);
    end
end
