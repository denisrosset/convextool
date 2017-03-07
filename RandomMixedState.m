function rho = RandomMixedState(dims)
% RandomMixedState([dA dB ...]) Generates a random mixed state of dimension dA x dB x ...
    d = prod(dims);
    rho = rand(d,d)*2-1 + 1i*(rand(d,d)*2-1);
    rho = rho*rho';
    rho = rho/trace(rho);
end
