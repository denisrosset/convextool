function rho = SchmidtNumberConeC(dims, k, def)
% SchmidtNumberConeC Approx. formulation of the cone of bipartite operators with bounded Schmidt number
%
% INPUTS
% dims      dims = [dA dB], dimensions of the A and B subsystems
% k         Upper bound on the operator Schmidt rank
% def       Definition of the separable cone with dimensions [dA*k dB*k]
%
% OUTPUTS
% rho       Complex (dA*dB)x(dA*dB) matrix representing the AB system
%           The A,B basis ordering is such that a product state
%           rho = rhoA (x) rhoB = kron(rhoA, rhoB)
    assert(def.dims(1) == dims(1)*k);
    assert(def.dims(2) == dims(2)*k);
    dA = dims(1);
    dB = dims(2);
    cvx_begin set sdp
    variable tau(dA*dB*k*k, dA*dB*k*k) hermitian
    variable rho(dA*dB, dA*dB)
    tau == SeparableConeC(def)
    tau1 = permute(reshape(tau, [dB k k dA dB k k dA]), [1 4 5 8 2 6 3 7]);
    I = eye(k*k);
    rho == reshape(reshape(tau1, dA*dB*dA*dB, k*k*k*k) * I(:), dA*dB, dA*dB);
    cvx_end
end
