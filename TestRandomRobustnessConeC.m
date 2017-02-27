for dA = 2:3
    for dB = 2:3
        if dA*dB < 9 % when the PPT criterion is exact
            pureState = (rand(dA*dB, 1)*2 - 1) + 1i * (rand(dA*dB, 1)*2 - 1);
            pureState = pureState / norm(pureState);
            exactValue = RandomRobustnessPureState(pureState, [dA dB]);                               
            % to compare to the conic optimization value
            density = pureState * pureState';
            cvx_solver sdpt3
            cvx_begin sdp
            variable nu nonnegative
            def = SymmetricExtensionDef([dA dB], 'exact');
            {nu density} == RandomRobustnessConeC(def);
            minimize nu
            cvx_end
            assert(abs(exactValue - nu) < 1e-7);
        end
    end
end
