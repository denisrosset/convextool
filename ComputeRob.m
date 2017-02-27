options = sdpsettings('verbose', 0, 'solver', 'sdpt3');
dA = 2;
dB = 3;
d = dA*dB;
ket = RandomPureState([dA dB]);
rho = ket*ket';
rho = (rho+rho')/2;
cvx_solver sdpt3
for k = 2:4
    %    inner = SymmetricExtensionDef([dA dB], 'inner', k, 'ppt', 'navascues', 'useSym', 1);
    outer = SymmetricExtensionDef([dA dB], 'outer', k, 'ppt', 'doherty', 'useSym', 1);
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
    variable subRhoS(d, d) hermitian
    minimize nu
    rho + nu*eye(d)/d == SymmetricExtensionConeC(outer)
    nu == trace(subRhoS)
    cvx_end
    lb(k) = nu;
end
