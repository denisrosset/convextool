function R = realify(C)
% REALIFY - convert a to be semipositive matrix from complex to real
    R = kron(real(C), eye(2)) + kron(imag(C), [0 -1; 1 0]);
end
