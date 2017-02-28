function Pabxy = MDIEW_Simulate(X, Y, rho, A, B)
    nX = length(X);
    nY = length(Y);
    nA = length(A);
    nB = length(B);
    Pabxy = zeros(nA, nB, nX, nY);
    for x = 1:length(X)
        for y = 1:length(Y)
            bigState = kron(kron(X{x}, rho), Y{y});
            for a = 1:length(A)
                for b = 1:length(B)
                    bigPOVM = kron(A{a}, B{b});
                    Pabxy(a,b,x,y) = real(trace(bigState*bigPOVM)); % remove imag. residuals
                end
            end
        end
    end
end
