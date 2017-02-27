function rrob = RandomRobustnessPureState(pureState, dims)
% RandomRobustnessPureState Computes the random robustness of a pure bipartite state
%
% The general definition of robustness is given just after Eq. (6)
% in http://link.aps.org/doi/10.1103/PhysRevA.59.141
% and the analytic formula for bipartite pure state in Eq.(41) of the same paper.
%
% INPUTS
% pureState      Complex vector of the coefficients of the pure state, ordered
%                such that a pure product state pureA (x) pureB would be written
%                pureState = kron(pureA, pureB)
% dims = [dA dB] Dimensions of the subsystems dA and dB
    dA = dims(1);
    dB = dims(2);
    % the Schmidt coefficients are given in decreasing order by svd in Matlab
    schmidt = svd(reshape(pureState, dB, dA));
    % analytic expression
    rrob = dA*dB*schmidt(1)*schmidt(2);
end
