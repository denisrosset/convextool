function S = SedumiOriginalSettings
% The default settings for SeDuMi provided by the interface YALMIP
% are taken from the default of an old SeDuMi version. Here are the
% defaults from the user manual of SeDuMi version 1.1.
    S = sdpsettings('solver', 'sedumi', 'sedumi.eps', 5e-8, ...
                    'sedumi.cg.qprec', 1, 'sedumi.cg.maxiter', 49, ...
                    'sedumi.stepdif', 2 , 'savesolveroutput', 1);
end
