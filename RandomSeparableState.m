function rho = RandomSeparableState(dims)
% RandomSeparableState([dA dB ...]) Generates a random density matrix corresponding to a separable multipartite state
    rho = zeros(prod(dims), prod(dims));
    for k = 1:prod(dims)
        contrib = 1;
        for i = 1:length(dims)
            contrib = kron(contrib, RandomMixedState(dims(i)));
        end
        rho = rho + contrib;
    end
    rho = rho / trace(rho);
end
