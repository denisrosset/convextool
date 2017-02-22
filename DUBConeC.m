function set = DUBConeC(dims)
% DUBConeC Entanglement cone of the SDP distillable entanglement upper bound
%
% {nu rho} = DUBConeC([dA dB]) returns the cone such that
%
% - rho is a bipartite density matrix, ordered such that the
%   product state rhoAB = rhoA (x) rhoB = kron(rhoA, rhoB)
% - nu is a nonnegative variable such that nu >= dub(rho)
%   where dub(rho) is given by the SDP of Eq. (8) of
%   http://link.aps.org/doi/10.1103/PhysRevA.94.050301
    dA = dims(1);
    dB = dims(2);
    d = dA*dB;
    cvx_begin set sdp
    variable rho(d, d) hermitian
    variable U(d, d) hermitian
    variable V(d, d) hermitian
    variable nu nonnegative
    M = reshape(U - V, [dB dA dB dA]);
    M = permute(M, [3 2 1 4]);
    M = reshape(M, [d d]);
    M - rho >= 0
    nu >= trace(U + V)
    U >= 0
    V >= 0
    rho >= 0
    cvx_end
    set = cvxtuple(struct('nu', nu, 'rho', rho));
end
