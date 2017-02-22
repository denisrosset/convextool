function set = NegativityConeC(dims)
% NegativityConeC Negativity entanglement measure cone
%
% {nu rho} = NegativityConeC([dA dB]) returns the cone such that
%
% - rho is a bipartite density matrix, ordered such that the
%   product state rhoAB = rhoA (x) rhoB = kron(rhoA, rhoB)
% - nu is a nonnegative variable such that nu >= negativity(rho)
%
% For rho not normalized, we interpret nu >= trace(rho)*negativity(rho/trace(rho))
    assert(length(dims) == 2);
    dA = dims(1);
    dB = dims(2);
    d = dA*dB;
    cvx_begin set sdp
    variable rho(d, d) hermitian
    variable sigmaP(d, d) hermitian
    variable sigmaM(d, d) hermitian
    variable nu nonnegative
    trace(sigmaM) <= nu
    rhoTA = reshape(rho, [dB dA dB dA]);
    rhoTA = permute(rhoTA, [3 2 1 4]);
    rhoTA = reshape(rhoTA, [d d]);
    sigmaP - sigmaM == rhoTA
    sigmaP >= 0
    sigmaM >= 0
    cvx_end
    set = cvxtuple(struct('nu', nu, 'rho', rho));
end
