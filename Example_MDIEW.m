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

% measurements
A = BellMeasurement;
B = BellMeasurement;

Pabxy = MDIEW_Simulate(X, Y, rho, A, B);

% serious stuff here
nX = length(X);
nY = length(Y);
nA = length(A);
nB = length(B);

% serious stuff begins here
cvx_clear
cvx_begin

variable Pi(dX*dY,dX*dY,nA,nB) hermitian % effective POVM elements
variable ent
variable nu(nA,nB) nonnegative
dual variable beta

minimize ent

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

beta: Pideal == Pabxy

sum_nu = 0;

for a = 1:nA
    for b = 1:nB
        {nu(a,b) Pi(:,:,a,b)} == NegativityConeC([dX dY])
        sum_nu = sum_nu + nu(a,b);
    end
end

ent == sum_nu

cvx_end