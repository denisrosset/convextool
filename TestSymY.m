options = sdpsettings('verbose', 0, 'solver', 'sdpt3');
rho0 = eye(9)/9;
R = (rand(9,9)*2-1) + 1i*(rand(9,9)*2-1);
rho1 = R*R';
rho1 = rho1/trace(rho1);
v = sdpvar;
rhov = rho0*(1-v) + rho1*v;
coeffsv = CoeffsFromOperator2(rhov, 3, 3);
values = [];
for useSym = 0:1
    for toReal = 0:1
        def = SymmetricExtensionDef([3 3], 2, 'ppt', 'doherty', 'useSym', useSym, 'toReal', toReal);
        Cons = SymmetricExtensionConeDualY(coeffsv, def);
        optimize(Cons, -v, options);
        values = [values double(v)];
    end
    def = SymmetricExtensionDef([3 3], 2, 'ppt', 'doherty', 'useSym', useSym, 'toReal', 0);
    Cons = SymmetricExtensionConePrimalY(rhov, def);
    optimize(Cons, -v, sdpsettings(options, 'dualize', 1));
    values = [values double(v)];
end
assert(max(values) - min(values) < 1e-5);
