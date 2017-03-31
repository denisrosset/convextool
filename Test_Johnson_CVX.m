% Tests the results in
% http://journals.aps.org/pra/pdf/10.1103/PhysRevA.88.032323
% Compatible quantum correlations: Extension problems for Werner and isotropic states
% Peter D. Johnson and Lorenza Viola
singlet = [0 0 0 0
           0 1 -1 0
           0 -1 1 0
           0 0 0 0]/2;
rho0 = eye(4)/4;
ks = [2 3 4];       % number of copies
vs = [2/3 5/9 1/2]; % visibilities

% use the bipartite code
for i = 1:3
    def = SeparableConeDef([2 2], 'outer', ks(i), 'ppt', []);
    cvx_clear
    cvx_begin sdp 
        variable v;
        maximize(v)
        subject to
        rho0*(1-v) + singlet*v == SeparableConeC(def)
    cvx_end
    assert(abs(v - vs(i)) < 1e-5);
end

% using the multipartite code in the bipartite case
for i = 1:3
    def = MultiSeparableConeDef([2 2], [1 ks(i)]);
    cvx_clear
    cvx_begin sdp
        variable v;
        maximize(v)
        subject to
        rho0*(1-v) + singlet*v == MultiSeparableConeC(def)
    cvx_end
    assert(abs(v - vs(i)) < 1e-5);
end
