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

settings = sdpsettings; %('removeequalities', 1);
% using the bipartite definition
def = SeparableConeDef([3 3], 'outer', 2, 'ppt', 'doherty');
yalmip('clear');
v = sdpvar;
CONS = SeparableConeY(def, 2*psiplus/7 + v * sigmaplus/7 + (5-v)*VsigmaplusV/7);
solvesdp(CONS, -v, settings);
v = double(v);
assert(abs(v - 3) < 1e-4);

yalmip('clear');
CONS = SeparableConeY(def, 2*psiplus/7 + v * sigmaplus/7 + (5-v)*VsigmaplusV/7);
solvesdp(CONS, v, settings);
v = double(v);
assert(abs(v - 2) < 1e-4);

% using the multipartite definition
def = MultiSeparableConeDef([3 3], [1 3], [0 1; 0 2; 0 3]);
yalmip('clear');
CONS = MultiSeparableConeY(def, 2*psiplus/7 + v * sigmaplus/7 + (5-v)*VsigmaplusV/7);
solvesdp(CONS, -v, settings);
v = double(v);
assert(abs(v - 3) < 1e-4);

yalmip('clear');
CONS = MultiSeparableConeY(def, 2*psiplus/7 + v * sigmaplus/7 + (5-v)*VsigmaplusV/7);
solvesdp(CONS, v, settings);
v = double(v);
assert(abs(v - 2) < 1e-4);
