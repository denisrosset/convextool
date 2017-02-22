function set = RandomRobustnessConeC(dims, k, ppt, useSym)
% RandomRobustnessConeC Random robustness entanglement measure cone
%
% {nu rho} = RandomRobustnessConeC([dA dB]) returns the cone such that
%
% - rho is a bipartite density matrix, ordered such that th
%   product state rhoAB = rhoA (x) rhoB = kron(rhoA, rhoB)
% - nu is a nonnegative variable such that nu >= randomrobustness(rho)
%   where randomrobustness(rho) is given by the Definition after Eq. (6)
%   of http://link.aps.org/doi/10.1103/PhysRevA.59.141
    assert(length(dims) == 2);
    dA = dims(1);
    dB = dims(2);
    d = dA*dB;
    rhoI = eye(d)/d;
    cvx_begin set sdp
    variable rho(d, d) hermitian
    variable nu nonnegative
    rho + nu * rhoI == SymmetricExtensionConeC(dims, k, ppt, useSym)
    cvx_end
    set = cvxtuple(struct('nu', nu, 'rho', rho));
end
