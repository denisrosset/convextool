classdef TestMeasuresApproximateCVX < CVXTestCase
% Tests the entanglement measures on the approximate case where
% when we use the Doherty hierarchy to approximate separable states
    properties (TestParameter)
        cone = {...
            {@AbsoluteRobustnessPureState, @AbsoluteRobustnessConeC} ...
            {@RandomRobustnessPureState, @RandomRobustnessConeC} ...
            {@GeneralizedRobustnessPureState, @GeneralizedRobustnessConeC} ...
            {@NegativityPureState, @NegativityConeC}
               };
        dA = {2 3 4};
        dB = {2 3 4};
    end
    methods(Test)
        function testNegativitySinglet(testCase)
            singlet = [0 0 0 0
                       0 1 -1 0
                       0 -1 1 0
                       0 0 0 0]/2;
        end
        function testPureState(testCase, cone, dA, dB)
            tol = ConfigTestOption('tolscale') * 1e-6;
            cbp = ConfigTestOption('cvx_begin_param');
            PureValue = cone{1};
            CvxCone = cone{2};
            def = SeparableConeDef([dA dB], 'outer', 2, 'ppt', 'doherty');
            pureState = (rand(dA*dB, 1)*2 - 1) + 1i * (rand(dA*dB, 1)*2 - 1);
            pureState = pureState / norm(pureState);
            exactValue = PureValue(pureState, [dA dB]);
            % to compare to the conic optimization value
            density = pureState * pureState';
            cvx_begin('sdp', cbp{:})
            variable nu nonnegative
            {nu density} == CvxCone([dA dB], def);
            minimize nu
            cvx_end
            assert(nu <= exactValue + tol);
        end
        function testSeparableState(testCase, cone, dA, dB)
            tol = ConfigTestOption('tolscale') * 1e-6;
            cbp = ConfigTestOption('cvx_begin_param');
            CvxCone = cone{2};
            def = SeparableConeDef([dA dB], 'outer', 2, 'ppt', 'doherty');
            density = RandomSeparableState([dA dB]);
            cvx_begin('sdp', cbp{:})
            variable nu nonnegative
            {nu density} == CvxCone([dA dB], def);
            minimize nu
            cvx_end
            assert(abs(nu) < tol);
        end
    end
end
