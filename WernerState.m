function rho = WernerState(v)
    singlet = [0  0  0 0
               0  1 -1 0
               0 -1  1 0
               0  0  0 0]/2;
    whiteNoise = eye(4)/4;
    rho = v*singlet + (1-v)*whiteNoise;
end
