options = sdpsettings('verbose', 0, 'solver', 'sdpt3');
rho0 = eye(4)/4;
R = (rand(4,4)*2-1) + 1i*(rand(4,4)*2-1);
rho1 = R*R';
rho1 = rho1/trace(rho1);
v = sdpvar;
rhov = rho0*(1-v) + rho1*v;
coeffsv = CoeffsFromOperator2(rhov, 2, 2);
values = [];
for realify = 0:1
    for useSym = 0:1
        Cons = SymmetricExtensionConeDual(coeffsv, 4, 'doherty', useSym, realify);
        optimize(Cons, -v, options);
        values = [values double(v)];
        Cons = SymmetricExtensionConePrimal(rhov, 2, 2, 4, 'doherty', useSym, realify);
        optimize(Cons, -v, sdpsettings(options, 'dualize', 1));
        values = [values double(v)];
    end
end
assert(max(values) - min(values) < 1e-5);
