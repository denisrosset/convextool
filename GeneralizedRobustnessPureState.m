function genrob = GeneralizedRobustnessPureState(pureState, dims)
% GeneralizedRobustnessPureState Computes the generalized robustness of a pure bipartite state
%
% Reduces to the robustness, see Theorem 3 of 
% http://journals.aps.org/pra/pdf/10.1103/PhysRevA.67.054305
    genrob = RobustnessPureState(pureState, dims);
end
