classdef TestSchmidtNumberConeCVX < CVXTestCase
    properties (TestParameter)
        dA = {3};
        dB = {3};
        kState = {1 2};
        kCone = {1 2};
    end
    methods(Test)
        function testSchmidtNumberCone(testCase, dA, dB, kState, kCone)
            if kState <= kCone
                tol = ConfigTestOption('tolscale') * 1e-6;
                cbp = ConfigTestOption('cvx_begin_param');
                ketAB = RandomSchmidtNumberPureState([dA dB], kState);
                rhoAB = ketAB*ketAB';
                rhoAB = (rhoAB + rhoAB')/2;
                def = SeparableConeDef([dA*kCone dB*kCone], 'outer', 2, 'ppt', []);
                cvx_begin('sdp', cbp{:})
                variable t
                minimize t 
                t*eye(dA*dB) + rhoAB == SchmidtNumberConeC([dA dB], kCone, def);
                cvx_end
                assert(t <= tol);
            end
        end
    end    
end
