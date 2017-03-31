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

% using the bipartite code
def = SeparableConeDef([2 2], 'outer', 2, 'ppt', []);
cvx_clear
cvx_begin sdp quiet
    variable v nonnegative
    maximize(v)
    subject to
    variable rhoAB(4,4) hermitian
    rhoAB >= 0
    rhoAB == state00*(1-v) + psiplus*v
    rhoAB == SeparableConeC(def);
cvx_end
assert(abs(v - 2/3) < 1e-5);

% using the multipartite code
def = MultiSeparableConeDef([2 2], [1 2], []);
cvx_clear
cvx_begin sdp quiet
    variable v nonnegative
    maximize(v)
    subject to
    variable rhoAB(4,4) hermitian
    rhoAB >= 0
    rhoAB == state00*(1-v) + psiplus*v
    rhoAB == MultiSeparableConeC(def);
cvx_end
assert(abs(v - 2/3) < 1e-5);
