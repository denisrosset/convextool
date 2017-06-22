classdef TestMeasuresExactCVX < CVXTestCase
% Tests the entanglement measures on the exact case where the PPT
% criterion is sufficient to detect separability
    properties (TestParameter)
        cone = {...
            {@AbsoluteRobustnessPureState, @AbsoluteRobustnessConeC} ...
            {@RandomRobustnessPureState, @RandomRobustnessConeC} ...
            {@GeneralizedRobustnessPureState, @GeneralizedRobustnessConeC} ...
            {@NegativityPureState, @NegativityConeC}
               };
        dA = {2 3};
        dB = {2 3};
    end
    methods(Test)
        function testSingletNegativity(testCase)
            tol = ConfigTestOption('tolscale') * 1e-6;
            cbp = ConfigTestOption('cvx_begin_param');
            singlet = [0 0 0 0; 0 1 -1 0; 0 -1 1 0; 0 0 0 0]/2;
            cvx_begin('sdp', cbp{:})
            variable nu nonnegative
            variable rho(4, 4) hermitian
            {nu rho} == NegativityConeC([2 2])
            rho == singlet
            minimize nu
            cvx_end
            assert(abs(nu - 1/2) < tol);
        end
        function testPureState(testCase, cone, dA, dB)
            tol = ConfigTestOption('tolscale') * 1e-6;
            cbp = ConfigTestOption('cvx_begin_param');
            PureValue = cone{1};
            CvxCone = cone{2};
            if dA*dB < 9 % when the PPT criterion is exact
                         % generate a random qubit-qubit pure state
                pureState = (rand(dA*dB, 1)*2 - 1) + 1i * (rand(dA*dB, 1)*2 - 1);
                pureState = pureState / norm(pureState);
                exactValue = PureValue(pureState, [dA dB]);
                % to compare to the conic optimization value
                density = pureState * pureState';
                cvx_begin('sdp', cbp{:})
                variable nu nonnegative
                {nu density} == CvxCone([dA dB]);
                minimize nu
                cvx_end
                assert(abs(exactValue - nu) < tol);
            end
        end
        function testSeparableState(testCase, cone, dA, dB)
            tol = ConfigTestOption('tolscale') * 1e-6;
            cbp = ConfigTestOption('cvx_begin_param');
            CvxCone = cone{2};
            if dA*dB < 9 % when the PPT criterion is exact
                         % generate a random qubit-qubit pure state
                density = RandomSeparableState([dA dB]);
                cvx_begin('sdp', cbp{:})
                variable nu nonnegative
                {nu density} == CvxCone([dA dB]);
                minimize nu
                cvx_end
                assert(abs(nu) < tol);
            end
        end
    end
end
