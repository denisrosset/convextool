singlet = [0 0 0 0
           0 1 -1 0
           0 -1 1 0
           0 0 0 0]/2;
cvx_solver sdpt3
cvx_begin sdp quiet
variable nu nonnegative
variable rho(4, 4) hermitian
{nu rho} == NegativityConeC([2 2])
rho == singlet
minimize nu
cvx_end
nu
