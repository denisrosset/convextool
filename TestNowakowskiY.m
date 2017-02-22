% Tests the example in
% Symmetric extendibility of quantum states
% Marcin  L.  Nowakowski
% http://iopscience.iop.org/article/10.1088/1751-8113/49/38/385301/meta;jsessionid=4DE611655AB4744ED75739E172E56C7C.ip-10-40-1-105
options = sdpsettings('verbose', 1, 'solver', 'sedumi');
state00 = [1 0 0 0
           0 0 0 0
           0 0 0 0
           0 0 0 0];
psiplus = [0 0 0 0
           0 1 1 0
           0 1 1 0
           0 0 0 0]/2;
v = sdpvar;
rhov = state00*(1-v) + psiplus*v;
coeffsv = CoeffsFromOperator2(rhov, 2, 2);
for useSym = 0:1
    for realify = 0:1
        Cons = SymmetricExtensionConeDualY(coeffsv, 2, [], useSym, realify);
        optimize(Cons, -v, options);
        assert(abs(double(v) - 2/3) < 1e-5);
    end
    Cons = SymmetricExtensionConePrimalY(rhov, [2 2], 2, [], useSym);
    optimize(Cons, -v, sdpsettings(options, 'dualize', 1));
    assert(abs(double(v) - 2/3) < 1e-5);
end
