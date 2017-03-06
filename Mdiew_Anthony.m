params = perfect_params;


inputs_order = [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 5 5 5 5 5 5 6 6 6 6 6 6
                1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 6 6 6 6 5 5 5 5 1 2 ...
                3 4 6 5 1 2 3 4 6 5 ];
Loss_correc = ((inputs_order(1,:)==5)+1).*((inputs_order(2,:)==5)+1).*((inputs_order(1,:)==6)+1).*((inputs_order(2,:)==6)+1);
Mat_loss_correc = transpose(repmat(Loss_correc,16,1));

i=1
data_file = dlmread(filename(i).name);
data_raw = data_file(:,[7:10,13:16,19:22,25:28]);
noise =[repmat((data_file(:,11)),1,4),repmat((data_file(:,17)),1,4),repmat((data_file(:,23)),1,4),repmat((data_file(:,29)),1,4)];
%correction of the loss due to Z basis measurement
clicks_i_ab = data_raw.*Mat_loss_correc - noise.*(Mat_loss_correc-1);
clicks_i_ab(clicks_i_ab<0)=0;
P_a_b_i = estimate_prob(clicks_i_ab, params);
[comp_a_b_i, sumsq] = make_prob_consistent(P_a_b_i, params);


for d=1:length(inputs_order)
    Pabxy_exp(:,:,inputs_order(1,d),inputs_order(2,d)) = comp_a_b_i(:,:,d)
end

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
                    Pideal(a,b,x,y) = trace(inputStates * Pi(:,:,a,b));
                end
            end
        end
    end

    % and here is the constraint: the effective POVM reproduces the correlations
    beta: Pideal == Pabxy_exp

    % now, let's find the lower bound on the entanglement measure
    sum_nu = 0;

    for a = 1:nA
        for b = 1:nB
            % constraint: entanglement cone
            {nu(a,b) Pi(:,:,a,b)} == NegativityConeC([dX dY])
            % replace by GeneralizedRobustnessConeC([dX dY])
            % or any of the other entanglement measure cones
            sum_nu = sum_nu + nu(a,b)/dX/dY;
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

dot(beta(:), Pabxy_exp(:))
