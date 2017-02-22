for dA = 2:3
    for dB = 2:3
        if dA*dB < 9 % when the PPT criterion is exact
            % generate a random qubit-qubit pure state
            pureState = (rand(dA*dB, 1)*2 - 1) + 1i * (rand(dA*dB, 1)*2 - 1);
            pureState = pureState / norm(pureState);
            exactValue = GeneralizedRobustnessPureState(pureState, [dA dB]);
            % to compare to the conic optimization value
            density = pureState * pureState';
            cvx_solver sdpt3
            cvx_begin sdp
            variable nu nonnegative
            def = SymmetricExtensionDef([dA dB], 1, 'ppt', 'doherty');
            {nu density} == GeneralizedRobustnessConeC(def);
            minimize nu
            cvx_end
            assert(abs(exactValue - nu) < 1e-8);
        end
    end
end
