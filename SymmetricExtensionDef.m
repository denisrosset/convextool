classdef SymmetricExtensionDef
    properties(SetAccess = immutable)
        dims;      % = [dA dB], dimensions of the subsystems composing the symmetric cone        

        dA;        % Dimension of the first subsystem A
        
        dB;        % Dimension of the second subsystem B
        
        approx;    % Either 'inner', 'outer' or 'approx'
                   % 'outer' is the Doherty-type symmetric extension
                   % 'exact' is valid only for dA*dB <= 6, and sets
                   %         k = 1 and 'ppt' = [1]
                   % 'inner' is the Navascues inner approximation, valid
                   %         when 'ppt' = []
                   %         or 'ppt' = 'navascues', that is [ceil(k/2)]

        k;         % Number of copies of subsystem B in the extension
        
        ppt;       % PPT constraints to add. Can be of the following form:
                   %  = []            no PPT constraints
                   %  = [t_1 ... t_n] with 1 <= t_j <= k (duplicates ignored)
                   %                  for each t_j, adds a PPT constraint such that a number t_j of copies
                   %                  of B is transposed
                   %  = 'doherty'     equivalent to [1 2 ... k], 
                   %                  PPT conditions in the original 2004 Doherty paper
                   %  = 'navascues'   equivalent to [ceil(k/2)], 
                   %                  PPT condition in the 2009 Navascues paper,
                   %                  also the way it is implemented in QETLAB, as of end 2016
                   % Default: [] (do not add PPT constraints)
        
        pptCuts;   % integer vector enumerating the cuts corresponding to the PPT constraints
                   % When ppt = 'doherty' or 'navascues', see interpretation above,
                   % otherwise pptCuts = ppt.
        
        useSym;    % Whether to use the symmetric subspace (default: 1)
    end
    methods
        function def = SymmetricExtensionDef(dims, approx, k, varargin)
        % SymmetricExtensionDef Defines a symmetric extension cone
        %
        % INPUTS
        %
        % dims       = [dA dB] Dimension of subsystems A and B
        % approx     = 'exact', 'outer', 'inner', see property definitions above
        % k          = number of copies of subsystem B in the extension
        % 
        % Additional properties are given each by a key/value pair. Possible keys
        % are 'ppt', 'useSym', and are explained in the 'properties' section of
        % 'SymmetricExtensionDef.m'
        %
        % Examples
        %
        % The symmetric extension considered in Doherty 2004, Section VII.A, consists
        % of the symmetric extension of a qutrit-qutrit state, using two copies of B, 
        % with two additional PPT constraints
        %
        % def = SymmetricExtensionDef([3 3], 2, 'outer', 'ppt', 'doherty')
        %
        % which is equivalent to
        %
        % def = SymmetricExtensionDef([3 3], 2, 'outer', 'ppt', [1 2])
        %
        % To use a formulation without symmetry handling, write:
        %   
        % def = SymmetricExtensionDef([3 3], 2, 'outer', 'ppt', 'doherty', 'useSym', 0)
        %
        % To obtain the PPT exact formulation for dA * dB <= 6, write:
        % def = SymmetricExtensionDef([dA dB], 1, 'exact')
        %
        % References
        % Navascues 2009, DOI: 10.1103/PhysRevLett.103.160404
        % Doherty 2004, DOI: 10.1103/PhysRevA.69.022308
            assert(length(dims) == 2);
            def.dims = dims(:)';
            def.dA = def.dims(1);
            def.dB = def.dims(2);
            switch approx
              case 'exact'
                if def.dA * def.dB > 6 && def.dA > 1 && def.dB > 1
                    error('Exact formulations of the symmetric cone are valid for 2x2 and 2x3 states');
                end
                if nargin < 3
                    k = 1;
                end
                def.k = k;
                def.approx = approx;
                def.ppt = [1];
                def.pptCuts = [1];
              otherwise
                def.k = k;
                def.approx = approx;
                def.ppt = [];
                def.pptCuts = [];
                def.useSym = 1;
            end
            assert(def.k >= 1);
            
            % Interpret varargin

            assert(mod(length(varargin), 2) == 0, 'Each key should have an associated value');

            toString = @(var) evalc(['disp(var)']);

            for keyIndex = 1:2:length(varargin)
                key = varargin{keyIndex};
                value = varargin{keyIndex + 1};
                switch key
                  case 'ppt'
                    def.ppt = value;
                    if isequal(value, 'doherty')
                        def.pptCuts = 1:def.k;
                    elseif isequal(value, 'navascues')
                        def.pptCuts = ceil(def.k/2);
                    else
                        def.pptCuts = unique(value(:)');
                        assert(all(def.pptCuts >= 1) && all(def.pptCuts <= def.k));
                    end
                  case 'useSym'
                    if ~isequal(value, 0) && ~isequal(value, 1)
                        error('Option useSym should be 0/false or 1/true');
                    end
                    def.useSym = value;
                  otherwise
                    warning(['Unsupported option: ' toString(varargin{1})]);
                end
            end

            % verifies the sanity of the definition
            switch def.approx
              case 'exact'
                if length(def.pptCuts) == 0
                    error('Exact formulations of the symmetric cone require PPT constraints.');
                end
              case 'inner'
                if length(def.pptCuts) == 0
                    if def.k == 1
                        error('The inner approximation requires either PPT constraints or k > 1.')
                    end
                elseif isequal(def.pptCuts, ceil(def.k/2))
                    % good
                else
                    error('The inner approximation is only valid for ''ppt'' = [] or ''navascues''');
                end
            end
        end
        
        function [ARho ATau] = ConstraintRepresentsSym(def)
        % equality constraints that express that
        % tauFull(dB^k*dA, dB^k*dA) == rho(dB*dA, dB*dA)
        % with tau represented in the symmetric subspace
            
            dA = def.dA;
            dB = def.dB;
            k = def.k;
            
            
        end
        
    end % methods
end % classdef
