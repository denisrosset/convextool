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

options = sdpsettings('verbose', 1, 'solver', 'sedumi');

for ext = 2:7
    def = SeparableConeDef([3 3], 'outer', ext, 'ppt', 'doherty');
    yalmip clear;
    v = sdpvar;
    rhov = 2*psiplus/7 + v * sigmaplus/7 + (5-v)*VsigmaplusV/7;
    Cons = SeparableConeY(def, rhov);
    [model,recoverymodel,diagnostic,internalmodel] = export(Cons, v, options);
    At = model.A;
    b = model.b;
    c = model.C;
    K = model.K;
    save(['yalmip_good_doherty_' num2str(ext) '.mat'], 'At', 'b','c', 'K');
end
