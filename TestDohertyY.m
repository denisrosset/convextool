proj = @(x) x(:)*x(:)';
psiplus = [1 0 0 0 1 0 0 0 1]';
psiplus = proj(psiplus)/3;
C1 = [0 1 0
      0 0 0
      0 0 0];
C2 = [0 0 0
      0 0 1
      0 0 0];
C3 = [0 0 0
      0 0 0
      1 0 0];
sigmaplus = proj(C1(:))/3 + proj(C2(:))/3 + proj(C3(:))/3;
C1 = C1';
C2 = C2';
C3 = C3';
v = sdpvar;
VsigmaplusV = proj(C1(:))/3 + proj(C2(:))/3 + proj(C3(:))/3;
rhov = 2*psiplus/7 + v * sigmaplus/7 + (5-v)*VsigmaplusV/7;
coeffsv = CoeffsFromOperator2(rhov, 3, 3);
options = sdpsettings('solver', 'sdpt3', 'verbose', 0);
for useSym = 0:1
    for toReal = 0:1
        % dual formulation
        def = SymmetricExtensionDef([3 3], 2, 'ppt', 'doherty', 'useSym', useSym, 'toReal', toReal);
        Cons = SymmetricExtensionConeDualY(coeffsv, def);
        optimize(Cons, -v, options);
        assert(abs(double(v) - 3) < 1e-4);
        optimize(Cons, v, options);
        assert(abs(double(v) - 2) < 1e-4);
    end
    % primal formulation
    def = SymmetricExtensionDef([3 3], 2, 'ppt', 'doherty', 'useSym', useSym, 'toReal', 0);
    Cons = SymmetricExtensionConePrimalY(rhov, def);
    optimize(Cons, -v, sdpsettings(options, 'dualize', 1));
    assert(abs(double(v) - 3) < 1e-5);
    optimize(Cons, v, sdpsettings(options, 'dualize', 1));
    assert(abs(double(v) - 2) < 1e-5);
end
