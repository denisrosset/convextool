% Implements the multipartite separability example of Doherty, Parillo, Spedalieri, 2005
% DOI 10.1103/PhysRevA.71.032333
ket0 = [1;0];
ket1 = [0;1];
ketp = [1;1]/sqrt(2);
ketm = [1;-1]/sqrt(2);

psi{1} = kron(ket0, kron(ket1, ketp));
psi{2} = kron(ket1, kron(ketp, ket0));
psi{3} = kron(ketp, kron(ket0, ket1));
psi{4} = kron(ketm, kron(ketm, ketm));

rho = eye(8, 8);
for i = 1:4
    rho = rho - psi{i} * psi{i}';
end
rho = rho / 4;

% all the ppt cuts
cuts = [1 0 0
        2 0 0
        0 1 0
        1 1 0
        2 1 0];
        
def = MultiSeparableConeDef([2 2 2], [2 1 1], cuts);
dtau = def.symmetricExtensionSize;
cvx_clear
cvx_begin sdp
    variable tau(dtau, dtau) hermitian
    variable nu nonnegative
    minimize nu
    rho + eye(8) * nu  == MultiSeparableConeC(def);
cvx_end
    
assert(nu > 0);
