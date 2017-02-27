options = sdpsettings('verbose', 0, 'solver', 'sdpt3');
dA = 3;
dB = 3;
d = dA*dB;
ket = RandomPureState([dA dB]);
rho = ket*ket';
rho = (rho+rho')/2;
noise = (1-2*rand(d,d)) + 1i*(1-2*rand(d,d));
noise = noise'*noise;
noise = noise / trace(noise);
rho = 0.9*rho + 0.2*noise;
cvx_solver scs
for k = 1:5
    %    inner = SymmetricExtensionDef([dA dB], 'inner', k, 'ppt', 'navascues', 'useSym', 1);
    outer = SymmetricExtensionDef([dA dB], 'outer', k, 'ppt', 'doherty');
    %    cvx_begin sdp
    %    variable nu nonnegative
    %    variable subRhoS(d, d) hermitian
    %    minimize nu
    %    rho + nu*eye(d)/d == SymmetricExtensionConeC(inner)
    %    nu == trace(subRhoS)
    %    cvx_end
    %    ub(k) = nu;
    cvx_begin sdp
    variable nu nonnegative
    %    variable subRhoS(d, d) hermitian
    minimize nu
    {nu rho} == AbsoluteRobustnessConeC(outer);
    %    rho + nu*eye(d)/d == SymmetricExtensionConeC(outer)
    %    nu == trace(subRhoS)
    cvx_end
    lb(k) = nu;
end
