function tests = TestSymmetricSubspace
    tests = functiontests(localfunctions);
end

function testOneCopy(testCase)
    n = 1;
    for d = randi(200, 100, 1)'
        ind = randi(d, 100, 1);
        sub = SymmetricSubspace(d, n).indToSubSym(ind);
        assert(isequal(sub, ind));
        ind = SymmetricSubspace(d, n).subToIndSym(sub);
        assert(isequal(sub, ind));
    end
end

function testDimensionOne(testCase)
    d = 1;
    for n = randi(200, 100, 1)'
        m = randi(100);
        ind = ones(m, 1);
        sub = SymmetricSubspace(d, n).indToSubSym(ind);
        assert(isequal(sub, ones(m, n)));
    end
end

function testIndToSubFull(testCase)
    for i = 1:10
        d = randi(8);
        n = randi(8);
        ind = randi(d^n, 200, 1);
        sub = SymmetricSubspace(d, n).indToSubFull(ind);
        ind1 = SymmetricSubspace(d, n).subToIndFull(sub);
        assert(isequal(ind, ind1));
    end
end

function testSubToIndFull(testCase)
    for i = 1:10
        d = randi(8);
        n = randi(8);
        sub = randi(d, 200, n);
        ind = SymmetricSubspace(d, n).subToIndFull(sub);
        sub1 = SymmetricSubspace(d, n).indToSubFull(ind);
        assert(isequal(sub, sub1));
    end
end

function testSubToIndSym(testCase)
    for i = 1:10
        d = randi(6);
        n = randi(6);
        ss = SymmetricSubspace(d, n);
        fullSub = ss.indToSubFull(1:d^n);
        symSub = unique(sort(fullSub')', 'rows');
        dim = size(symSub, 1);
        assert(isequal(dim, ss.dim));
        test = randi(dim, 100, 1);
        subTest = symSub(test, :);
        symInd = ss.subToIndSym(subTest);
        assert(isequal(symInd, test));
    end
end

function testIndToSubSym(testCase)
    for i = 1:10
        d = randi(6);
        n = randi(6);
        ss = SymmetricSubspace(d, n);
        fullList = ss.indToSubFull(1:d^n);
        % list of canonical subindices, sorted in increasing order for each row
        canonical = unique(sort(fullList')', 'rows');
        sizeBasis = size(canonical, 1);
        assert(isequal(sizeBasis, ss.dim));
        test = randi(sizeBasis, 100, 1); % canonical indices to test
        testSub = ss.indToSubSym(test);
        correctSub = canonical(test, :);
        assert(isequal(testSub, correctSub));
    end
end

function testDimension(testCase)
    for i = 1:10
        d = randi(6);
        n = randi(6);
        ss = SymmetricSubspace(d, n);
        fullList = ss.indToSubFull(1:d^n);
        % list of canonical subindices, sorted in increasing order for each row
        canonical = unique(sort(fullList')', 'rows');
        sizeBasis = size(canonical, 1);
        assert(isequal(sizeBasis, ss.dim));
    end
end
