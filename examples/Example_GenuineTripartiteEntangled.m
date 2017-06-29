% unnormalized GHZ state
un_ghz = [1 0 0 0 0 0 0 1]';
rho = un_ghz * un_ghz'/2;
rho0 = eye(8)/8;
def = SeparableConeDef([4 2], 'outer', 2, 'ppt', 'doherty');
cvx_solver mosek
cvx_begin sdp
variable t nonnegative
variable sepAB_C(8,8) hermitian semidefinite
variable sepAC_B(8,8) hermitian semidefinite
variable sepBC_A(8,8) hermitian semidefinite
sepAB_C == SeparableConeC(def)
sepAC_B == SeparableConeC(def)
sepBC_A == SeparableConeC(def)  
%                                 C A B C A B     C B A C B A
comp1 = permute(reshape(sepAB_C, [2 2 2 2 2 2]), [1 3 2 4 6 5]);
%                                 B A C B A C    
comp2 = permute(reshape(sepAC_B, [2 2 2 2 2 2]), [3 1 2 6 4 5]);
%                                 A B C A B C
comp3 = permute(reshape(sepBC_A, [2 2 2 2 2 2]), [3 2 1 6 5 4]);

comp1 = reshape(comp1, [8 8]);
comp2 = reshape(comp2, [8 8]);
comp3 = reshape(comp3, [8 8]);
t * rho + (1-t)*rho0 == comp1 + comp2 + comp3
maximize t
cvx_end
