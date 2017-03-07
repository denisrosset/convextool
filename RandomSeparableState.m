function rhoAB = RandomSeparableState(dims)
% RandomSeparableState([dA dB]) Generates a random density matrix corresponding to a separable bipartite state
    dA = dims(1);
    dB = dims(2);
    rhoAB = zeros(dA*dB, dA*dB);
    for k = 1:dA*dB
        rhoAB = rhoAB + rand * kron(RandomMixedState(dA), ...
                                    RandomMixedState(dB));
    end
    rhoAB = rhoAB / trace(rhoAB);
end
