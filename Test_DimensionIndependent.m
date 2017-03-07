% Tests that the generalized and absolute robustness are
% independent of the dimension
dA = 2;
dB = 2;
pureState = RandomPureState([dA dB]);
% to compare to the conic optimization value
rhoAB = pureState * pureState';
sigmaAB = RandomSeparableState([dA dB]);
def = SeparableConeDef([dA^2 dB^2], 'outer', 2, 'ppt', 'doherty');
tauAB = reshape(kron(rhoAB, sigmaAB), [dA dB dA dB dA dB dA dB]);
tauAB = permute(tauAB, [1 3 2 4 5 7 6 8]);
tauAB = reshape(tauAB, [dA^2*dB^2 dA^2*dB^2]);

cvx_solver sdpt3

exactValue = GeneralizedRobustnessPureState(pureState, [dA dB]);
cvx_begin sdp
variable nu nonnegative
{nu tauAB} == GeneralizedRobustnessConeC([dA^2 dB^2], def);
minimize nu
cvx_end
assert(abs(exactValue - nu) < 1e-5);

exactValue = AbsoluteRobustnessPureState(pureState, [dA dB]);
cvx_begin sdp
variable nu nonnegative
{nu tauAB} == AbsoluteRobustnessConeC([dA^2 dB^2], def);
minimize nu
cvx_end
assert(abs(exactValue - nu) < 1e-5);
