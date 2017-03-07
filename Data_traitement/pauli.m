function M = pauli(k)
% PAULI - Pauli matrices
%
% k = 1,2,3 for x,y,z
    switch k
      case 1
        M = [0 1; 1 0];
      case 2
        M = [0 -1i; 1i 0];
      case 3
        M = [1 0; 0 -1];
      otherwise
        M = [1 0; 0 1];
    end
end
