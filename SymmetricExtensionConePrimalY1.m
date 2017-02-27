function [Cons Info] = SymmetricExtensionConePrimalY1(rhoAB, def)
% SymmetricExtensionConePrimalY Compute SDP constraints corresponding to symmetric extension cones
%
% Formulation in the YALMIP primal canonical form (= using equalities)
% To gain advantage of this formulation, set sdpsettings('dualize', 1)
%
% INPUTS
%
% rhoAB        Complex (dA*dB)x(dA*dB) matrix representing the AB system
%              The A,B basis ordering is such that a product state
%              rhoAB = rhoA (x) rhoB = kron(rhoA, rhoB)
% 
% def          Symmetric cone definition, see SymmetricExtensionDef
%
% OUTPUTS
%
% Cons         All the SDP constraints associated with the cone
%
% Info         TBD

%    if isequal(def.approx, 'inner')
%        error('Only implemented for the primal/CVX formulation SymmetricExtensionConeC.');
%    end
%    if ~def.useSym
%        warning('The nonsymmetric primal form is inefficient. Set ''useSym'', 1');
%    end
%    if def.toReal
%        warning('Real formulation is not supported in primals. Use SymmetricExtensionConeDualY.');
%    end
    dims = def.dims;
    k = def.k;
    ppt = def.ppt;
    useSym = def.useSym;
    dA = dims(1);
    dB = dims(2);
    d = dA * dB;
    assert(size(rhoAB, 1) == dA * dB);
    assert(size(rhoAB, 2) == dA * dB);
    usePPT = length(ppt) > 0;
    if def.toReal
        matrixType = 'symmetric';
        fieldDim = 2;
    else
        matrixType = 'hermitian';
        fieldDim = 1;
    end
    if useSym
        L = nchoosek(dB+k-1, k); % Eq. (25) of Doherty2004
        tauS = sdpvar(L*fieldDim, L*fieldDim, matrixType);
    else
        % constraint that the symmetric extension represents
        % the bipartite state
        reprN = d*(d+1)/2;
        reprI = 1;
        reprLHS = sparse(reprN, dA^2*dB^2);
        reprRHS = sparse(reprN, dA^2*dB^(2*k));
        for r = 1:dA*dB
            for c = r:dA*dB
                [rB rA] = ind2sub([dB dA], r);
                [cB cA] = ind2sub([dB dA], c);
                reprLHS(reprI, r + (c-1)*d) = 1;
                sz1 = [dB dB^(k-1) dA dB dB^(k-1) dA];
                for jr = 1:dB^(k-1)
                    reprRHS(ind, vsub2ind(sz1, [rB jr rA cB jr cA])) = 1;
                end
            end
        end
        % constraint that the symmetric extension is indeed symmetric
        factB = cumprod(dB*ones(1, k));
        factB = [1 factB(1:k-1)];
        s2iB = @(x) factB*(x(:)-1)+1; % sub2ind for B subsystems
        canIndB = SymmetricCanonicalIndices(dB, k);
        % this gives 
        
        error('Not implemented');
    end
end
