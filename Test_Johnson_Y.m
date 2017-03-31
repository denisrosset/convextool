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

settings = sdpsettings;

% use the bipartite code
for i = 1:3
    yalmip('clear');
    def = SeparableConeDef([2 2], 'outer', ks(i), 'ppt', []);
    v = sdpvar;
    CONS = SeparableConeY(def, rho0*(1-v) + singlet*v);
    solvesdp(CONS, -v, settings); % maximize v
    v = double(v);
    assert(abs(v - vs(i)) < 1e-5);
end

% using the multipartite code in the bipartite case
for i = 1:3
    yalmip('clear');
    def = MultiSeparableConeDef([2 2], [1 ks(i)]);
    v = sdpvar;
    CONS = MultiSeparableConeY(def, rho0*(1-v) + singlet*v);
    solvesdp(CONS, -v, settings); % maximize v
    v = double(v);
    assert(abs(v - vs(i)) < 1e-5);
end
