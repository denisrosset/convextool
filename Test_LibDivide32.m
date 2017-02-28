function tests = TestLibDivide
    tests = functiontests(localfunctions);
end

function testCorrect(testCase)
    IMAX = (2^32) - 1;
    for D = randi(IMAX, 1, 100)
        ld = LibDivide32(D);
        n = randi(IMAX, 1, 100);
        res1 = ld.div(n);
        res2 = floor(n/D);
        assert(isequal(res1, res2));
    end
end
