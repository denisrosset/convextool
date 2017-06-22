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

inputsX = [pX(:) mX(:) pY(:) mY(:) pZ(:) mZ(:)];
inputsX = reshape(inputsX, [2 2 6]);
inputsY = inputsX;

B = BellMeasurement;
nX = 6;
nY = nX;
inputsX = inputsX(:,:,1:nX);
inputsY = inputsY(:,:,1:nY);
nB = 4;
PQ = zeros(nB,nX,nY);
for b = 1:nB
    for x = 1:nX
        for y = 1:nY
            PQ(b,x,y) = trace(kron(inputsX(:,:,x), inputsY(:,:,y))*B{b});
        end
    end
end
P0 = ones(nB,nX,nY)/nB;
cvx_clear
cvx_begin sdp
variable t
maximize t
(t*PQ + (1-t)*P0) == TemporalMDIEWSetC(nB, inputsX, inputsY)
cvx_end
t