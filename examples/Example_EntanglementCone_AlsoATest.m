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

% formulating and solving the problem
cvx_clear
cvx_begin sdp

    variable nu nonnegative
    variable rho(4,4) hermitian
    dual variable W
    
    minimize nu

    subject to
    W: rho == werner
    {nu rho} == RandomRobustnessConeC([dA dB]);
    % try also
    % {nu rho} == AbsoluteRobustnessConeC([dA dB]);
    % {nu rho} == GeneralizedRobustnessConeC([dA dB]);
    % {nu rho} == DUBConeC([dA dB]);
    % {nu rho} == NegativityConeC([dA dB]);
    
cvx_end

disp('The Werner state with parameter')
v

disp('has random robustness')
nu

disp('with quantitative entanglement witness')
W
