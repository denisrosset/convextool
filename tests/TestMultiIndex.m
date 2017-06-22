function tests = TestMultiIndex
    tests = functiontests(localfunctions);
end

function testOneDim(testCase)
    n = 1;
    for d = randi(200, 1, 100)
        dims = [d];
        ind = randi(d, 100, 1);
        sub = MultiIndex(dims).indToSub(ind);
        assert(isequal(sub, ind));
        ind = MultiIndex(dims).subToInd(sub);
        assert(isequal(sub, ind));
    end
end

function testDimensionOne(testCase)
    for n = randi(200, 1, 100)
        dims = ones(1, n);
        m = randi(100);
        ind = ones(m, 1);
        sub = MultiIndex(dims).indToSub(ind);
        assert(isequal(sub, ones(m, n)));
    end
end

function TestZeroDim(testCase)
    dims = [];
    m = randi(100);
    ind = ones(m, 1);
    sub = MultiIndex(dims).indToSub(ind);
    assert(isequal(size(sub), [m 0]));
    ind1 = MultiIndex(dims).subToInd(sub);
    assert(isequal(ind, ind1));
end

function TestSubToInd(testCase)
    for n = randi([2 8], 1, 100) % sub2ind does not support dims=[d1]
        dims = randi(8, 1, n);
        m = randi(100);
        csub = cell(1, n); % cell array of subindices
        msub = ones(m, n); % matrix of subindices
        for c = 1:n
            col = randi(dims(c), m, 1);
            msub(:, c) = col;
            csub{c} = col;
        end
        ind1 = MultiIndex(dims).subToInd(msub);
        ind2 = sub2ind(dims, csub{:});
        assert(isequal(ind1, ind2));
    end
end

function TestIndToSub(testCase)
    for n = randi([2 8], 1, 100)
        dims = randi(8, 1, n);
        m = randi(100);
        ind = randi(prod(dims), m, 1);
        csub = cell(1, n);
        [csub{:}] = ind2sub(dims, ind);
        msub = MultiIndex(dims).indToSub(ind);
        for c = 1:n
            assert(isequal(csub{c}, msub(:,c)));
        end
    end
end
