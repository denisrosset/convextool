% generate a random qubit-qubit pure state
purestate = (rand(4, 1)*2 - 1) + 1i * (rand(4, 1)*2 - 1);
purestate = purestate / norm(purestate);

% find the Schmidt coefficients
schmidt = svd(reshape(purestate, 2, 2));
n1 = 2;
n2 = 2;
a1 = schmidt(1);
a2 = schmidt(2);
randomrobustness = n1*n2*a1*a2; % eq.41 of 
density = purestate * purestate';
cvx_solver sdpt3
cvx_begin sdp
variable nu nonnegative
variable rho(4, 4) hermitian
{nu density} == RandomRobustnessConeC([2 2], 1, 1, 'doherty')
minimize nu
cvx_end
assert(abs(randomrobustness - nu) < 1e-8);
