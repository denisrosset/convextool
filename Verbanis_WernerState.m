function rho = WernerState(v)
    singlet = [1  0  0 -j
               0  0 0 0
               0 0  0 0
               j 0  0 1]/2;
    whiteNoise = eye(4)/4;
    rho = v*singlet + (1-v)*whiteNoise;
end
