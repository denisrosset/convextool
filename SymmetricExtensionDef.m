function def = SymmetricExtensionDef(dims, k, varargin)
% SymmetricExtensionDef Defines a symmetric extension cone
%
% INPUTS
%
% dims       = [dA dB]  Dimension of subsystems A and B
% k          Number of copies of subsystem B in the extension
%
% followed by the key/value pairs:
%
% 'ppt'      PPT constraints to add. Can be of the following form:
%            = []            no PPT constraints
%            = [t_1 ... t_n] with 1 <= t_j <= k (duplicates ignored)
%                            for each t_j, adds a PPT constraint such that a number t_j of copies
%                            of B is transposed
%            = 'doherty'     equivalent to [1 2 ... k], 
%                            PPT conditions in the original 2004 Doherty paper
%            = 'navascues'   equivalent to [ceil(k/2)], 
%                            PPT condition in the 2009 Navascues paper,
%                            also the way it is implemented in QETLAB, as of end 2016
%            Default: [] (do not add PPT constraints)
%
% 'useSym'   Whether to use the symmetric subspace (default: 1)
%
% 'toReal'   Whether to use a real SDP formulation (default: 0), however
%            real formulations are only supported in YALMIP's dual formulations.
%            For the primal YALMIP&CVX formulations this option is ignored.
%
% Examples
%
% The symmetric extension considered in Doherty 2004, Section VII.A, consists
% of the symmetric extension of a qutrit-qutrit state, using two copies of B, 
% with two additional PPT constraints
%
% def = SymmetricExtensionDef([3 3], 2, 'ppt', 'doherty')
%
% which is equivalent to
%
% def = SymmetricExtensionDef([3 3], 2, 'ppt', [1 2])
%
% To use a formulation without symmetry handling, write:
%   
% def = SymmetricExtensionDef([3 3], 2, 'ppt', 'doherty', 'useSym', 0)
%
% To enable real formulations:
% def = SymmetricExtensionDef([3 3], 2, 'ppt', 'doherty', 'toReal', 1)
%
%
% References
% Navascues 2009, DOI: 10.1103/PhysRevLett.103.160404
% Doherty 2004, DOI: 10.1103/PhysRevA.69.022308
    assert(length(dims) == 2);
    assert(k >= 1);
    defaults = struct('dims', dims(:)', ...
                 'k', k, ...
                 'ppt', [], ...
                 'useSym', 1, ...
                 'toReal', 1);
    def = interpret(defaults, varargin{:});
    function def = interpret(def, varargin)
        toString = @(var) evalc(['disp(var)']);
        switch length(varargin)
          case 0
          case 1
            error(['Missing value for option' toString(varargin{1})]);
          otherwise
            key = varargin{1};
            value = varargin{2};
            switch key
              case 'ppt'
                if isequal(value, 'doherty')
                    def.ppt = 1:def.k;
                elseif isequal(value, 'navascues')
                    def.ppt = ceil(def.k/2);
                else
                    def.ppt = unique(value(:)');
                    assert(all(def.ppt >= 1) && all(def.ppt <= def.k));
                end
              case 'useSym'
                if ~isequal(value, 0) && ~isequal(value, 1)
                    error('Option useSym should be 0/false or 1/true');
                end
                def.useSym = value;
              case 'toReal'
                if ~isequal(value, 0) && ~isequal(value, 1)
                    error('Option toReal should be 0/false or 1/true');
                end
                def.toReal = value;
              otherwise
                warning(['Unsupported option: ' toString(varargin{1})]);
            end
            def = interpret(def, varargin{3:end});
        end
    end
end
