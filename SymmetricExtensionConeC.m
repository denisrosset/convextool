function rhoAB = SymmetricExtensionConeC(def)
% SymmetricExtensionConeC Approximation of the cone of separable operators
%
% Formulation in the CVX standard form (= SeDuMi primal, using equalities)
% There is no support of "realification". 
%
% INPUTS
% def        Symmetric cone definition, see SymmetricExtensionDef
%
% OUTPUTS
% rhoAB      Complex (dA*dB)x(dA*dB) matrix representing the AB system
%            The A,B basis ordering is such that a product state
%            rhoAB = rhoA (x) rhoB = kron(rhoA, rhoB)
    
    if ~def.useSym
        warning('The nonsymmetric primal form is inefficient. Set ''useSym'', 1')
    end
    if def.toReal
        warning('Real formulation is not supported in primals. Use SymmetricExtensionConeDualY.');
    end
    dims = def.dims;
    k = def.k;
    ppt = def.ppt;
    usePPT = length(ppt) > 0;
    useSym = def.useSym;
    dA = dims(1);
    dB = dims(2);
    
    if k == 1 % handle separately the PPT criterion alone, without extension
        d = dA*dB;
        if isequal(ppt, [])
            warning('The symmetric 1-extension without PPT constraint is trivial');
            cvx_begin set sdp
            variable rhoAB(d, d) semidefinite hermitian % main variable
            rhoAB >= 0
            cvx_end
        else
            cvx_begin set sdp
            variable rhoAB(d, d) semidefinite hermitian % main variable
            variable pptAB(d, d) semidefinite hermitian % partial transpose
            rhoTA = reshape(rhoAB, [dB dA dB dA]);
            rhoTA = permute(rhoTA, [3 2 1 4]);
            rhoTA = reshape(rhoTA, [d d]);
            rhoTA == pptAB
            rhoAB >= 0
            pptAB >= 0
            cvx_end
        end
        return
    end

    % start the convex cone definition
    cvx_begin set sdp
    variable rhoAB(dA*dB, dA*dB) hermitian % main variable

    if useSym
        [~, G] = BasisSymmetricSubspace(dB, k);
        dBext = size(G, 2);
        variable tau(dA*dBext, dA*dBext) hermitian % variable
        tau >= 0;
        conv = kron(eye(dA), G);
        tauFull = conv * tau * conv';
    else
        dBext = dB^k;
        [~, nPi, dPi] = ProjectorSymmetricSubspace(dB, k);
        nPi = kron(eye(dA), nPi);
        variable tauFull(dA*dBext, dA*dBext) hermitian % variable
        tauFull >= 0;
        nPi * tauFull * nPi == dPi * dPi * tauFull; % symmetry
    end
    tauAB = reshape(tauFull, [dB dB^(k-1) dA dB dB^(k-1) dA]);
    tauAB = permute(tauAB, [1 3 4 6 2 5]);
    tauAB = reshape(tauAB, [dB*dA*dB*dA dB^(k-1)*dB^(k-1)]);
    % partial trace
    tauAB = tauAB * reshape(eye(dB^(k-1)), dB^(k-1)*dB^(k-1), 1);
    tauAB = reshape(tauAB, [dB*dA dB*dA]);
    tauAB == rhoAB; % partial trace reproduces rhoAB
    if usePPT
        if useSym
            for i = 1:length(ppt)
                k1 = ppt(i);
                k2 = k - ppt(i);
                [~, G1] = BasisSymmetricSubspace(dB, k1);
                [~, G2] = BasisSymmetricSubspace(dB, k2);
                dB1 = size(G1, 2);
                dB2 = size(G2, 2);
                conv = kron(eye(dA), kron(G1, G2) \ G); % split the symmetric subspace
                tauPPT = conv * tau * conv';
                tauPPT = reshape(tauPPT, [dB1 dB2 dA dB1 dB2 dA]);
                tauPPT = permute(tauPPT, [4 2 3 1 5 6]);
                tauPPT = reshape(tauPPT, [dB1*dB2*dA dB1*dB2*dA]);
                variable(['tauPPTvar' num2str(i) '(dB1*dB2*dA, dB1*dB2*dA)'], 'hermitian');
                tauPPTvar = eval(['tauPPTvar' num2str(i)]);
                tauPPTvar == tauPPT; % PPT constraint
                tauPPTvar >= 0
            end
        else
            for i = 1:length(ppt)
                k1 = ppt(i);
                k2 = k - ppt(i);
                tauPPT = reshape(tauFull, [dB^k1 dB^k2 dA dB^k1 dB^k2 dA]);
                tauPPT = permute(tauPPT, [4 2 3 1 5 6]);
                tauPPT = reshape(tauPPT, [dB^k*dA dB^k*dA]);
                variable(['tauPPTvar' num2str(i) '(dB^k*dA, dB^k*dA)'], 'hermitian');
                tauPPTvar = eval(['tauPPTvar' num2str(i)]);
                tauPPTvar == tauPPT; % PPT constraint
                tauPPTvar >= 0
            end
        end
    end
    cvx_end
end
