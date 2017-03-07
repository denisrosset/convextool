singlet = [0 0 0 0
           0 1 -1 0
           0 -1 1 0
           0 0 0 0]/2;
cvx_clear
cvx_solver sdpt3
cvx_begin sdp quiet
variable nu nonnegative
variable rho(4, 4) hermitian
{nu rho} == NegativityConeC([2 2])
rho == singlet
minimize nu
cvx_end
assert(abs(nu - 1/2) < 1e-7);

density = RandomSeparableState([2 2]);
cvx_solver sdpt3
cvx_begin sdp quiet
    variable nu nonnegative
    {nu density} == NegativityConeC([dA dB]);
    minimize nu
cvx_end
assert(nu < 1e-8);
