function set = GeneralizedRobustnessConeC(def)
% GeneralizedRobustnessConeC Outer approximation of the generalized robustness entanglement measure cone
%
% {nu rho} = GeneralizedRobustnessConeC([dA dB]) returns the cone such that
%
% - rho is a bipartite density matrix, ordered such that the
%   product state rhoAB = rhoA (x) rhoB = kron(rhoA, rhoB)
% - nu is a nonnegative variable such that nu >= genrob(rho)
%   where genrob(rho) is given by replacing rhoS \in SepCone by
%   rhoS >= 0 in Eq. (11) of http://link.aps.org/doi/10.1103/PhysRevA.59.141,
%   or, alternatively, following http://journals.aps.org/pra/pdf/10.1103/PhysRevA.67.054305
%
% The robustness is computed with respect to an approximation of
% the separable cone given in the parameter 'def', obtained by SymmetricExtensionDef.
%
% For dA*dB <= 6, the approximation is exact if def.ppt = 'doherty' or [1].

    dims = def.dims;
    dA = dims(1);
    dB = dims(2);
    d = dA*dB;
    rhoI = eye(d)/d;
    cvx_begin set sdp
    variable rho(d, d) hermitian
    variable subRhoM(d, d) hermitian % subnormalized state rhoS
    variable nu nonnegative
    subRhoM >= 0
    rho + subRhoM == SymmetricExtensionConeC(def)
    nu == trace(subRhoM)
    cvx_end
    set = cvxtuple(struct('nu', nu, 'rho', rho));
end
