function rho = RandomSeparableState(dims, n)
% RandomSeparableState([dA dB ...]) Generates a random density
% matrix corresponding to a separable multipartite state
    if nargin < 2
        n = prod(dims);
    end
    rho = zeros(prod(dims), prod(dims));
    for k = 1:n
        contrib = 1;
        for i = 1:length(dims)
            contrib = kron(contrib, RandomPureState(dims(i)));
        end
        rho = rho + rand * contrib*contrib';
    end
    rho = rho / trace(rho);
    rho = (rho + rho')/2;
end
