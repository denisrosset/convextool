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
VsigmaplusV = proj(C1(:))/3 + proj(C2(:))/3 + proj(C3(:))/3;
rhov = 2*psiplus/7 + v * sigmaplus/7 + (5-v)*VsigmaplusV/7;
coeffsv = CoeffsFromOperator2(rhov, 3, 3);
options = sdpsettings(S, 'verbose', 1, 'solver', 'sdpt3');
Cons = SymmetricExtensionCone(coeffsv, 2, false, true);
optimize(Cons, -v, options);
double(v)
diag = optimize(Cons, v, options);
double(v)
