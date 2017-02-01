% Tests the results in
% http://journals.aps.org/pra/pdf/10.1103/PhysRevA.88.032323
% Compatible quantum correlations: Extension problems for Werner and isotropic states
% Peter D. Johnson and Lorenza Viola
options = sdpsettings('verbose', 0);
singlet = [0 0 0 0
           0 1 -1 0
           0 -1 1 0
           0 0 0 0]/2;
rho0 = eye(4)/4;
v = sdpvar;
rhov = rho0*(1-v) + singlet*v;
coeffsv = CoeffsFromOperator2(rhov, 2, 2);
ks = [2 3 4];       % number of copies
vs = [2/3 5/9 1/2]; % visibilities
for i = 1:3
    Cons = SymmetricExtensionCone(coeffsv, ks(i), false, false);
    optimize(Cons, -v, options);
    assert(abs(double(v) - vs(i)) < 1e-5);
    Cons = SymmetricExtensionCone(coeffsv, ks(i), true, false);
    optimize(Cons, -v, options);
    assert(abs(double(v) - vs(i)) < 1e-5);
end
