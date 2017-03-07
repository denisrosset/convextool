function [Hcomplex, Hreal, HcomplexPT, HrealPT] = hermitian2(x)
    d = 2;
    Hcomplex = zeros(d*d, d*d);
    Hreal = zeros(d*d*2, d*d*2);
    HcomplexPT = zeros(d*d, d*d);
    HrealPT = zeros(d*d*2, d*d*2);
    assert(length(x) == 16);
    for j = 1:4
        for k = 1:4
            l = sub2ind([4 4], j, k);
            Hcomplex = Hcomplex + x(l) * kron(pauli(j), pauli(k));
            Hreal = Hreal + x(l) * realify(kron(pauli(j), pauli(k)));
            HcomplexPT = HcomplexPT + x(l) * kron(pauli(j).', pauli(k));
            HrealPT = HrealPT + x(l) * realify(kron(pauli(j).', pauli(k)));
        end
    end
end
