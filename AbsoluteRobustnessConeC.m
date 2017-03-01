function set = AbsoluteRobustnessConeC(dims, def)
% AbsoluteRobustnessConeC Outer approximation of the absolute robustness entanglement measure cone
%
% {nu rho} = AbsoluteRobustnessConeC([dA dB], def) returns the cone such that
%
% - rho is a bipartite density matrix, ordered such that the
%   product state rhoAB = rhoA (x) rhoB = kron(rhoA, rhoB)
% - nu is a nonnegative variable such that nu >= robustness(rho)
%   where robustness(rho) is given by Eq. (11) of http://link.aps.org/doi/10.1103/PhysRevA.59.141
%
% The absolute robustness is computed with respect to a formulation of
% the separable cone given in the parameter 'def', obtained by SeparableConeDef.
%
% If 'def' is omitted, then it gets the default value SeparableConeDef(dims, 'exact'),
% which will only works for dA*dB <= 6.
    if nargin < 2
        def = SeparableConeDef(dims, 'exact');
    end
    assert(isequal(dims, def.dims));
    dA = dims(1);
    dB = dims(2);
    d = dA*dB;
    rhoI = eye(d)/d;
    cvx_begin set sdp
    variable rho(d, d) hermitian
    variable subRhoS(d, d) hermitian % subnormalized state rhoS
    variable nu nonnegative
    subRhoS == SeparableConeC(def)
    rho + subRhoS == SeparableConeC(def)
    nu == trace(subRhoS)
    cvx_end
    set = cvxtuple(struct('nu', nu, 'rho', rho));
end
