function neg = NegativityPureState(pureState, dims)
% NegativityPureState Computes the negativity of a pure bipartite state
%
% The general definition of negativity is given in Eq. (1) of
% http://dx.doi.org/10%2E1103/PhysRevA%2E65%2E032314
%
% INPUTS
% pureState      Complex vector of the coefficients of the pure state, ordered
%                such that a pure product state pureA (x) pureB would be written
%                pureState = kron(pureA, pureB)
% dims = [dA dB] Dimensions of the subsystems dA and dB
    assert(length(dims) == 2);
    dA = dims(1);
    dB = dims(2);
    % Schmidt coefficients
    s = svd(reshape(pureState, dB, dA));
    neg = (sum(s)^2 - 1)/2;
end
