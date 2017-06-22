function P = TemporalMDIEWSetC(nB, inputsX, inputsY, def)
    assert(length(size(inputsX)) == 3);
    assert(length(size(inputsY)) == 3);
    dX = size(inputsX, 1);
    assert(dX == size(inputsX, 2));
    nX = size(inputsX, 3);
    dY = size(inputsY, 1);
    assert(dY == size(inputsY, 2));
    nY = size(inputsY, 3);
    if nargin < 4
        def = SeparableConeDef([dX dY], 'exact');
    end
    dims = def.dims;
    dA = dims(1);
    dB = dims(2);
    d = dA*dB;
    cvx_begin set sdp
    variable P(nB, nX, nY)
    variable mu(dX*dY, dX*dY, nB) hermitian
    for b = 1:nB
        mu(:,:,b) == SeparableConeC(def)
        for x = 1:nX
            for y = 1:nY
                P(b,x,y) == trace(kron(inputsX(:,:,x), inputsY(:,:,y)) * mu(:,:,b))
            end
        end
    end
    cvx_end
end
