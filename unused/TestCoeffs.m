function tests = TestCoeffs
    tests = functiontests(localfunctions);
end

function testSinglet(testCase)
    singlet = [0 0 0 0
               0 1 -1 0
               0 -1 1 0
               0 0 0 0]/2;
    coeffs = CoeffsFromOperator2(singlet, 2, 2);
    H = OperatorFromCoeffs2(coeffs);
    assert(isequal(singlet, H));
end

function testNoise(testCase)
    noise = eye(4)/4;
    coeffs = CoeffsFromOperator2(noise, 2, 2);
    H = OperatorFromCoeffs2(coeffs);
    assert(isequal(noise, H));
end
