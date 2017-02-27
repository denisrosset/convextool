function sub = vind2sub(sz, ind)
    ind = ind - 1;
    d = length(sz);
    sub = zeros(1, length(sz));
    for i = 1:d
        sub(i) = mod(ind, sz(i));
        ind = ind - sub(i);
        sub(i) = sub(i) + 1;
        ind = ind / sz(i);
    end
end
