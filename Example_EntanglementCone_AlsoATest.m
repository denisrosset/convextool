% dimensions
dA = 2;
dB = 2;

% state
rho = [0  0  0  0
       0  1 -1  0
       0 -1  1  0
       0  0  0  0]/2;
rho0 = eye(4)/4;
v = 1/2;
werner = v * rho + (1 - v) * rho0;

% the random robustness uses a formulation of the separable cone,
% which we create here. It also provides the subsystem dimensions.
def = SeparableConeDef([dA dB], 'exact');

% formulating and solving the problem
cvx_clear
cvx_begin sdp

    variable nu nonnegative
    variable rho(4,4) hermitian
    dual variable W
    
    minimize nu

    subject to
    W: rho == werner
    {nu rho} == RandomRobustnessConeC(def);
    % try also
    % {nu rho} == AbsoluteRobustnessConeC(def);
    % {nu rho} == GeneralizedRobustnessConeC(def);
    % {nu rho} == DUBConeC([dA dB]);
    % {nu rho} == NegativityConeC([dA dB]);
    
cvx_end

disp('The Werner state with parameter')
v

disp('has random robustness')
nu

disp('with quantitative entanglement witness')
W
