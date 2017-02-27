for d = 1:5
    for k = 1:5
        [~, B1] = BasisSymmetricSubspace(d, k);
        [~, B2] = BasisSymmetricSubspace(d, k, 1);
        assert(norm(abs(B1 - B2), 'inf') == 0);
    end
end
