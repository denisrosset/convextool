s_x = [0 1; 1 0];
s_y = [0 -1i; 1i 0];
s_z = [1 0; 0 -1];
I2 = eye(2);
pX = (I2 + s_x)/2;
mX = (I2 - s_x)/2;
pY = (I2 + s_y)/2;
mY = (I2 - s_y)/2;
pZ = (I2 + s_z)/2;
mZ = (I2 - s_z)/2;

% Alice inputs
X = {pX mX pY mY pZ mZ};
% Bob inputs
Y = X;

% and their dimension
dX = 2;
dY = 2;

% underlying state
v = 1;
rho = WernerState(v);

% Alice and Bob perform a Bell measurement
% the variable A is a cell(1,4) with A{a} 4x4 Hermitian POVM elements
A = BellMeasurement;
B = BellMeasurement;

% simulates the correlations
Pabxy = MDIEW_Simulate(X, Y, rho, A, B);

nX = length(X);
nY = length(Y);
nA = length(A);
nB = length(B);

% serious stuff begins here
cvx_clear
cvx_begin

    variable Pi(dX*dY,dX*dY,nA,nB) hermitian % effective POVM elements
    variable entanglement_lb
    variable nu(nA,nB) nonnegative

    % dual variable = MDIEW
    dual variable beta

    minimize entanglement_lb

    subject to

    % we first compute the left hand side of the constraint
    expressions Pideal(nA,nB,nX,nY) nonnegative
    for x = 1:nX
        for y = 1:nY
            inputStates = kron(X{x}, Y{y});
            for a = 1:nA
                for b = 1:nB
                    Pideal(a,b,x,y) = dX*dY*trace(inputStates * Pi(:,:,a,b));
                end
            end
        end
    end

    % and here is the constraint: the effective POVM reproduces the correlations
    beta: Pideal == Pabxy

    % now, let's find the lower bound on the entanglement measure
    sum_nu = 0;

    for a = 1:nA
        for b = 1:nB
            % constraint: entanglement cone
            def = SeparableConeDef([dX dY], 'exact');
            {nu(a,b) Pi(:,:,a,b)} == NegativityConeC([dX dY])
            % replace by GeneralizedRobustnessConeC(SeparableConeDef([dX dY], 'exact'))
            % or any of the other entanglement measure cones
            sum_nu = sum_nu + nu(a,b);
        end
    end

    entanglement_lb == sum_nu

cvx_end

disp('A lower bound on the entanglement from the MDIEW correlations is:')

entanglement_lb

disp('')

disp('Whereas the entanglement in the state is')

cvx_clear
cvx_begin sdp quiet
variable statenu nonnegative
minimize statenu
{statenu rho} == NegativityConeC([dX dY])
cvx_end

statenu

disp('')

disp('MDIEW coefficients in dual variable beta')

disp('')

disp('Recomputing the MDIEW value')

dot(beta(:), Pabxy(:))