function rob = RobustnessPureState(pureState, dims)
% RobustnessPureState Computes the robustness of a pure bipartite state
%
% The general definition of robustness is given in 
% Eq. (11) of http://link.aps.org/doi/10.1103/PhysRevA.59.141
% and the analytic formula for bipartite pure state in Eq.(30) of the same paper.
%
% INPUTS
% pureState      Complex vector of the coefficients of the pure state, ordered
%                such that a pure product state pureA (x) pureB would be written
%                pureState = kron(pureA, pureB)
% dims = [dA dB] Dimensions of the subsystems dA and dB
% find the Schmidt coefficients
    dA = dims(1);
    dB = dims(2);
    schmidt = svd(reshape(pureState, dB, dA));
    rob = sum(schmidt)^2 - 1;
end
