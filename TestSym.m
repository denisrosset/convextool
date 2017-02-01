options = sdpsettings('verbose', 0, 'solver', 'sdpt3');
rho0 = eye(4)/4;
R = (rand(4,4)*2-1) + 1i*(rand(4,4)*2-1);
rho1 = R*R';
rho1 = rho1/trace(rho1);
v = sdpvar;
rhov = rho0*(1-v) + rho1*v;
coeffsv = CoeffsFromOperator2(rhov, 2, 2);
disp('Without symmetry basis')
Cons = SymmetricExtensionCone(coeffsv, 4, false, true);
optimize(Cons, -v, options);
vnosym = double(v);
disp('With symmetry basis')
Cons = SymmetricExtensionCone(coeffsv, 4, true, true);
optimize(Cons, -v, options);
vsym = double(v);
assert(abs(vnosym - vsym) < 1e-5);
