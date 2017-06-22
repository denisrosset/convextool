classdef TestDUBConeCCVX < CVXTestCase
    properties (TestParameter)
        dA = {2 3 4};
        dB = {2 3 4};
        a = {1/3 1/2 2/3};
    end
    methods(Test)
        function testExactValues(testCase, a)
        % Example of proposition 5 of http://link.aps.org/doi/10.1103/PhysRevA.94.050301
        % Here, we have EW = log2(nu);
        % Bob in rows, Alice in columns
            tol = ConfigTestOption('tolscale') * 1e-6;
            cbp = ConfigTestOption('cvx_begin_param');
            psi0 = zeros(3, 3);
            psi1 = zeros(3, 3);
            psi2 = zeros(3, 3);
            psi0(1,2) = sqrt(a);
            psi0(2,1) = sqrt(1-a);
            psi1(1,3) = sqrt(a);
            psi1(3,1) = sqrt(1-a);
            psi2(2,3) = sqrt(a);
            psi2(3,2) = sqrt(1-a);
            psi0 = psi0(:);
            psi1 = psi1(:);
            psi2 = psi2(:);
            rho = (psi0*psi0' + psi1*psi1' + psi2*psi2')/3;
            exactValue = 1 + sqrt(a*(1-a));
            cvx_begin('sdp', cbp{:})
            variable nu nonnegative
            {nu rho} == DUBConeC([3 3]);
            minimize nu
            cvx_end
            assert(abs(nu - exactValue) < tol);
        end
        function testSeparableStates(testCase, dA, dB)
            tol = ConfigTestOption('tolscale') * 1e-6;
            cbp = ConfigTestOption('cvx_begin_param');
            density = RandomSeparableState([dA dB]);
            cvx_begin('sdp', cbp{:})
            variable nu nonnegative
            {nu density} == DUBConeC([dA dB]);
            minimize nu
            cvx_end
            assert(log2(nu) < tol);
        end
    end    
end
