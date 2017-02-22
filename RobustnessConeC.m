function set = RobustnessConeC(def)
% RobustnessConeC Lower bound on the random robustness entanglement measure cone
%
% {nu rho} = RobustnessConeC([dA dB]) returns the cone such that
%
% - rho is a bipartite density matrix, ordered such that the
%   product state rhoAB = rhoA (x) rhoB = kron(rhoA, rhoB)
% - nu is a nonnegative variable such that nu >= robustness(rho)
%   where robustness(rho) is given by Eq. (11) of http://link.aps.org/doi/10.1103/PhysRevA.59.141
%
% The robustness is computed with respect to an approximation of
% the separable cone given in the parameter 'def', obtained by SymmetricExtensionDef.
    dims = def.dims;
    dA = dims(1);
    dB = dims(2);
    d = dA*dB;
    rhoI = eye(d)/d;
    cvx_begin set sdp
    variable rho(d, d) hermitian
    variable subRhoS(d, d) hermitian % subnormalized state rhoS
    variable nu nonnegative
    subRhoS == SymmetricExtensionConeC(def)
    rho + subRhoS == SymmetricExtensionConeC(def)
    nu == trace(subRhoS)
    cvx_end
    set = cvxtuple(struct('nu', nu, 'rho', rho));
end
