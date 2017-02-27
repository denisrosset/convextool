function ind = vsub2ind(sz, sub)
    d = length(sz);
    cum = 1;
    ind = 1;
    for i = 1:d
        ind = ind + cum * (sub(i) - 1);
        cum = cum * sz(i);
    end
end