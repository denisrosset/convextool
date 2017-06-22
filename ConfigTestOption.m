function value = ConfigTestOption(key)
    solver = 'mosek';
    verbose = 1;
    
    switch solver
      otherwise
        tolscale = 1;
    end
   
    switch key
      case 'cvx_begin_param'
        if verbose
            value = {};
        else
            value = {'quiet'};
        end
      case 'solver'
        value = solver;
      case 'verbose'
        value = 1;
      case 'tolscale'
        value = tolscale;
      otherwise
        error(['Unknown test option ' key]);
    end
end
