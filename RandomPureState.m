function ket = RandomPureState(dims)
    d = prod(dims);
    ket = (2*rand(d, 1) - 1) + 1i * (2*rand(d, 1) - 1);
    ket = ket/norm(ket);
end
