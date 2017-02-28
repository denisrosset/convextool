% Example of proposition 5 of http://link.aps.org/doi/10.1103/PhysRevA.94.050301
% Here, we have EW = log2(nu);
% Bob in rows, Alice in columns
for a = [1/3 1/2 2/3]
    a = 1/2;
    psi0 = zeros(3, 3);
    psi1 = zeros(3, 3);
    psi2 = zeros(3, 3);
    psi0(1,2) = sqrt(a);
    psi0(2,1) = sqrt(1-a);
    psi1(1,3) = sqrt(a);
    psi1(3,1) = sqrt(1-a);
    psi2(2,3) = sqrt(a);
    psi2(3,2) = sqrt(1-a);
    psi0 = psi0(:);
    psi1 = psi1(:);
    psi2 = psi2(:);
    rho = (psi0*psi0' + psi1*psi1' + psi2*psi2')/3;
    exactValue = 1 + sqrt(a*(1-a));
    cvx_clear
    cvx_begin sdp quiet
    variable nu nonnegative
    {nu rho} == DUBConeC([3 3]);
    minimize nu
    cvx_end
    assert(abs(nu - exactValue) < 1e-7);
end
