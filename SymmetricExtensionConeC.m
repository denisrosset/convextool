function rho = SymmetricExtensionConeC(def)
% SymmetricExtensionConeC Approximation of the cone of separable operators
%
% Formulation in the CVX standard form (= SeDuMi primal, using equalities)
% There is no support of "realification". 
%
% INPUTS
% def        Symmetric cone definition, see SymmetricExtensionDef
%
% OUTPUTS
% rho        Complex (dA*dB)x(dA*dB) matrix representing the AB system
%            The A,B basis ordering is such that a product state
%            rho = rhoA (x) rhoB = kron(rhoA, rhoB)
    
    if ~def.useSym
        error('For now, we only support the symmetric form.');
    end
    
    dims = def.dims;
    dA = dims(1);
    dB = dims(2);
    d = dA*dB;
    k = def.k;
    
    if k == 1 % handle separately the PPT criterion alone, without extension
        d = dA*dB;
        if length(ppt) == 0
            warning('The symmetric 1-extension without PPT constraint is trivial');
            cvx_begin set sdp
            variable rhoAB(d, d) semidefinite hermitian % main variable
            cvx_end
        else
            error('TODO: better implementation');
            %            cvx_begin set sdp
            %            variable rhoAB(d, d) semidefinite hermitian % main variable
            %            variable pptAB(d, d) semidefinite hermitian % partial transpose
            %            hermMatIndices = SymmetricSubspace(dA*dB, 2);
            %            hermMatIndices 
            %            rhoTA = reshape(rhoAB, [dB dA dB dA]);
            %            rhoTA = permute(rhoTA, [3 2 1 4]);
            %            rhoTA = reshape(rhoTA, [d d]);
            %            rhoTA == pptAB
            %            cvx_end
        end
        return
    end

    ppt = def.ppt;
    nPPT = length(ppt);
    tauS = SymmetricSubspace(dB, k);
    dBsym = tauS.dim;
    pptS1 = cell(nPPT, 1);
    pptS2 = cell(nPPT, 1);
    for i = 1:nPPT
        p = ppt(i);
        pptS1{i} = SymmetricSubspace(dB, p);
        pptS2{i} = SymmetricSubspace(dB, k - p);
    end
    
    % Constraints
    % tau reproduces rho
    
    % (full) index subspaces AB and R = B^(k-1)
    B_A = MultiIndex([dB dA]);
    AB_AB = MultiIndex([d d]);
    B_A_B_A = MultiIndex([dB dA dB dA]);
    Bk = MultiIndex(dB*ones(1, k));
    B_R = MultiIndex([dB dB^(k-1)]);
    Bsym_A = MultiIndex([dBsym dA]);
    Bsym_A_Bsym_A = MultiIndex([dBsym dA dBsym dA]);
    %     = MultiIndex([dB dB*ones(1, k-1)]);
    
    nRest = dB^(k-1); % dimension over which we perform the partial trace
    traceOver = (1:nRest)';
    nRepr = (d+1)*d/2; % upper triangle dimension
    reprRho = sparse(nRepr, d^2);
    reprTau = sparse(nRepr, dBsym*dA*dBsym*dA);
    i = 1;
    for r = 1:d
        rBA = B_A.indToSub(r);
        rB = rBA(:,1);
        rA = rBA(:,2);
        rBfull = Bk.indToSub(B_R.subToInd([rB*ones(nRest, 1) traceOver])); % format B_R to B_B.._B
        rBsym = tauS.subToIndSym(rBfull);
        for c = r:d % restrict constraints to upper triangle
            cBA = B_A.indToSub(c);
            cB = cBA(:,1);
            cA = cBA(:,2);
            cBfull = Bk.indToSub(B_R.subToInd([cB*ones(nRest, 1) traceOver])); % format B_R to B_B.._B
            cBsym = tauS.subToIndSym(cBfull);
            
            reprRho(i, AB_AB.subToInd([r c])) = 1;
            tauInd = Bsym_A_Bsym_A.subToInd([rBsym rA*ones(nRest, 1) cBsym cA*ones(nRest, 1)]);
            for j = 1:nRest
               reprTau(i, tauInd(j)) = reprTau(i, tauInd(j)) + 1;
            end
            i = i + 1;
        end
    end

    % start the convex cone definition
    
    cvx_begin set sdp
    variable rho(dB*dA, dB*dA) hermitian % main variable
    variable tau(dBsym*dA, dBsym*dA) hermitian % symmetric extension in symmetric basis
    rho >= 0
    tau >= 0
    reprRho * rho(:) == reprTau * tau(:)
    cvx_end
end
