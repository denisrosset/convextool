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
cvx_clear
cvx_solver mosek
for k = 5
    inner = SeparableConeDef([dA dB], 'inner', k, 'ppt', 'navascues');
    outer = SeparableConeDef([dA dB], 'outer', k, 'ppt', 'doherty');
    %    cvx_begin sdp
    %    cvx_precision low
    %    variable nu nonnegative
    %    minimize nu
    %    {nu rho} == RandomRobustnessConeC(inner);
    %    cvx_end
    %    ub(k) = nu;
    cvx_begin sdp
    cvx_precision low
    variable nu nonnegative
    minimize nu
    {nu rho} == RandomRobustnessConeC(outer);
    cvx_end
    lb(k) = nu;
end
