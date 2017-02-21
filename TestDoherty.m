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
options = sdpsettings('verbose', 1, 'solver', 'sedumi');
for useSym = 0:1
    for realify = 0:1
        % dual formulation
        Cons = SymmetricExtensionConeDual(coeffsv, 2, 'doherty', useSym, realify);
        optimize(Cons, -v, options);
        abs(double(v)-3)
        assert(abs(double(v) - 3) < 1e-4);
        diag = optimize(Cons, v, options);
        assert(abs(double(v) - 2) < 1e-4);
        abs(double(v)-2)
        % primal formulation
        Cons = SymmetricExtensionConePrimal(rhov, 3, 3, 2, 'doherty', useSym, realify);
        optimize(Cons, -v, sdpsettings(options, 'dualize', 1));
        abs(double(v)-3)
        assert(abs(double(v) - 3) < 1e-4);
        diag = optimize(Cons, v, sdpsettings(options, 'dualize', 1));
        assert(abs(double(v) - 2) < 1e-4);
        abs(double(v)-2)
    end
end
