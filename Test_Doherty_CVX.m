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
for useSym = 1%0:1
    def = SeparableConeDef([3 3], 'outer', 3, 'ppt', 'doherty', 'useSym', useSym);
    cvx_clear
    cvx_begin sdp quiet
    variable v
    maximize(v)
    subject to
    2*psiplus/7 + v * sigmaplus/7 + (5-v)*VsigmaplusV/7 == SeparableConeC(def);
    cvx_end
    assert(abs(v - 3) < 1e-4);
    
    cvx_clear
    cvx_begin sdp quiet
    variable v
    minimize(v)
    subject to
    2*psiplus/7 + v * sigmaplus/7 + (5-v)*VsigmaplusV/7 == SeparableConeC(def);
    cvx_end
    assert(abs(v - 2) < 1e-4);
end
