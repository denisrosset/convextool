function [Real Imag] = RealifiedParts(M)
    dim = size(M, 1)/2;
    M = reshape(M, [2 dim 2 dim]);
    M = permute(M, [2 4 1 3]);
    M = reshape(M, [dim*dim 4]);
    Real = M * reshape(eye(2), 4, 1)/2;
    Imag = M * reshape([0 -1; 1 0], 4, 1)/2;
    Real = reshape(Real, dim, dim);
    Imag = reshape(Imag, dim, dim);
end
