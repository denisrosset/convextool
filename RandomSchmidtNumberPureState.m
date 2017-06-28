function ket = RandomSchmidtNumberPureState(dims, k)
    assert(length(size(dims)) == 2);
    dA = dims(1);
    dB = dims(2);
    RA = (2*rand(dA,dA)-1) + 1i*(2*rand(dA,dA)-1);    
    [UA, DA] = eig(RA);
    RB = (2*rand(dB,dB)-1) + 1i*(2*rand(dB,dB)-1);    
    [UB, DB] = eig(RB);
    UA = UA(:,1:k);
    UB = UB(:,1:k);
    D = diag(rand(k, 1));
    ket = (UA*D*UB')';
    ket = ket(:);
    ket = ket/norm(ket);
end
