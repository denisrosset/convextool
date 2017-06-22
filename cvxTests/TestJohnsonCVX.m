classdef TestJohnsonCVX < CVXTestCase
    properties(TestParameter)
        k = {2 3 4}; % number of copies
    end
    methods(Test)
        function test(testCase, k)
        % Tests the results in
        % http://journals.aps.org/pra/pdf/10.1103/PhysRevA.88.032323
        % Compatible quantum correlations: Extension problems for Werner and isotropic states
        % Peter D. Johnson and Lorenza Viola
            tol = ConfigTestOption('tolscale') * 1e-6;
            cbp = ConfigTestOption('cvx_begin_param');
            vs = [2/3 5/9 1/2];
            realv = vs(k - 1); % visibilities
            singlet = [0 0 0 0
                       0 1 -1 0
                       0 -1 1 0
                       0 0 0 0]/2;
            rho0 = eye(4)/4;
            % use the bipartite code
            def = SeparableConeDef([2 2], 'outer', k, 'ppt', []);
            cvx_begin('sdp', cbp{:})
            variable v;
            maximize(v)
            subject to
            rho0*(1-v) + singlet*v == SeparableConeC(def)
            cvx_end
            assert(abs(v - realv) < tol);
            % using the multipartite code in the bipartite case
            def = MultiSeparableConeDef([2 2], [1 k]);
            cvx_begin('sdp', cbp{:})
            variable v;
            maximize(v)
            subject to
            rho0*(1-v) + singlet*v == MultiSeparableConeC(def)
            cvx_end
            assert(abs(v - realv) < tol);
        end
    end
end
