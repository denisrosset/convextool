function dm = qubit_state(vec)
% QUBIT_STATE Qubit density matrix from Bloch vector
    dm = (eye(2) + vec(1) * pauli(1) + vec(2) * pauli(2) + ...
         vec(3) * pauli(3))/2;
end
