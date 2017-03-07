% Tests that random separable states are inside the Doherty approximated cone
for dA = 2:5
    for dB = 2:5
        rhoAB = RandomSeparableState([dA dB]);
        def = SeparableConeDef([dA dB], 'outer', 2, 'ppt', 'doherty');
        cvx_solver mosek
        cvx_clear
        cvx_begin sdp
        variable t
        minimize t 
        t*eye(dA*dB) + rhoAB == SeparableConeC(def);
        cvx_end
        tol = 1e-6;
        assert(t <= tol);
    end
end
