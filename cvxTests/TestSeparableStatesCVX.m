classdef TestSeparableStatesCVX < CVXTestCase
    properties (TestParameter)
        dA = {2 3 4};
        dB = {2 3 4};
        
    end
    methods(Test)
        function testSeparableStates(testCase, dA, dB)
            tol = ConfigTestOption('tolscale') * 1e-6;
            cbp = ConfigTestOption('cvx_begin_param');
            rhoAB = RandomSeparableState([dA dB]);
            def = SeparableConeDef([dA dB], 'outer', 2, 'ppt', 'doherty');
            cvx_begin('sdp', cbp{:})
            variable t
            minimize t 
            t*eye(dA*dB) + rhoAB == SeparableConeC(def);
            cvx_end
            assert(t <= tol);
        end
    end    
end
