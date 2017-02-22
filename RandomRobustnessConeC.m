function set = RandomRobustnessConeC(def)
% RandomRobustnessConeC Outer approximation of the random robustness entanglement measure cone
%
% {nu rho} = RandomRobustnessConeC([dA dB]) returns the cone such that
%
% - rho is a bipartite density matrix, ordered such that the
%   product state rhoAB = rhoA (x) rhoB = kron(rhoA, rhoB)
% - nu is a nonnegative variable such that nu >= randomrobustness(rho)
%   where randomrobustness(rho) is given by the Definition after Eq. (6)
%   of http://link.aps.org/doi/10.1103/PhysRevA.59.141
%
% The random robustness is computed with respect to an approximation of the
% separable cone given in the parameter 'def', obtained by SymmetricExtensionDef.
%
% For dA*dB <= 6, the approximation is exact if def.ppt = 'doherty' or [1].
    dims = def.dims;
    dA = dims(1);
    dB = dims(2);
    d = dA*dB;
    rhoI = eye(d)/d;
    cvx_begin set sdp
    variable rho(d, d) hermitian
    variable nu nonnegative
    rho + nu * rhoI == SymmetricExtensionConeC(def)
    cvx_end
    set = cvxtuple(struct('nu', nu, 'rho', rho));
end
