classdef TestMultiSeparableStatesCVX < CVXTestCase
    properties (TestParameter)
        dims = {[2 2 3] [3 2 2] [2 3 2]};
    end
    methods(Test)
        function testSeparableStates(testCase, dims)
            tol = ConfigTestOption('tolscale') * 1e-6;
            cbp = ConfigTestOption('cvx_begin_param');
            rho = RandomSeparableState(dims, 1);
            cuts = [1 0 0
                    1 1 0
                    0 1 0];
            def = MultiSeparableConeDef(dims, [2 2 2], cuts);
            cvx_begin('sdp', cbp{:})
            variable t
            minimize t 
            t*eye(prod(dims)) + rho == MultiSeparableConeC(def);
            cvx_end
            assert(t <= tol);
        end
    end    
end
