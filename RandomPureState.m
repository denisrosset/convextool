function ket = RandomPureState(dims)
    dA = dims(1);
    dB = dims(2);
    d = dA*dB;
    ket = (2*rand(d, 1) - 1) + 1i * (2*rand(d, 1) - 1);
    ket = ket/norm(ket);
end
