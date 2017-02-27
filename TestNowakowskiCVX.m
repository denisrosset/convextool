% Tests the example in
% Symmetric extendibility of quantum states
% Marcin  L.  Nowakowski
% http://iopscience.iop.org/article/10.1088/1751-8113/49/38/385301/meta;jsessionid=4DE611655AB4744ED75739E172E56C7C.ip-10-40-1-105
state00 = [1 0 0 0
           0 0 0 0
           0 0 0 0
           0 0 0 0];
psiplus = [0 0 0 0
           0 1 1 0
           0 1 1 0
           0 0 0 0]/2;
for useSym = 1
    def = SymmetricExtensionDef([2 2], 'outer', 2, 'ppt', [], 'useSym', useSym);
    cvx_solver sdpt3
    cvx_begin sdp
    %cvx_dualize off
    variable v nonnegative
    maximize(v)
    subject to
    variable rhoAB(4,4) hermitian
    rhoAB >= 0
    rhoAB == state00*(1-v) + psiplus*v
    rhoAB == SymmetricExtensionConeC(def);
    cvx_end
    assert(abs(v - 2/3) < 1e-5);
end
