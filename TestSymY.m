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
    for realify = 0:1
        Cons = SymmetricExtensionConeDualY(coeffsv, 2, 'doherty', useSym, realify);
        optimize(Cons, -v, options);
        values = [values double(v)];
    end
    Cons = SymmetricExtensionConePrimalY(rhov, [3 3], 2, 'doherty', useSym);
    optimize(Cons, -v, sdpsettings(options, 'dualize', 1));
    values = [values double(v)];
end
assert(max(values) - min(values) < 1e-5);
