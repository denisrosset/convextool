function [Cons Info] = SymmetricExtensionConePrimalY(rhoAB, def)
% SymmetricExtensionConePrimalY Compute SDP constraints corresponding to symmetric extension cones
%
% Formulation in the YALMIP primal canonical form (= using equalities)
% To gain advantage of this formulation, set sdpsettings('dualize', 1)
% There is no support of "realification".
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
% Info         Extra information to use in diagnostics
%
% Info.tauC    SDP constraint of the symmetric extension tau >= 0
%
%              For 'useSym' = 0
% Info.symC    Symmetric equality constraint
% Info.tau =
% Info.tauFull Density matrix of the symmetric extension
%
%              For 'useSym' = 1
% Info.tau     Density matrix of the symmetric extension (reduced sym. basis)
% Info.tauFull Density matrix of the symmetric extension (after sym. basis expansion)
%
%              For length(def.ppt) > 1
% Info.tauPT   Cell array of matrices of partial transposes
% Info.tauPTC  Cell array of PPT constraints PT{i} >= 0
%
%              For k > 1
% Info.reprC   The constraint that trace_{B2..Bk}(tau) == rhoAB
     
    if ~def.useSym
        warning('The nonsymmetric primal form is inefficient. Set ''useSym'', 1');
    end
    if def.toReal
        warning('Real formulation is not supported in primals. Use SymmetricExtensionConeDualY.');
    end
    dims = def.dims;
    k = def.k;
    ppt = def.ppt;
    useSym = def.useSym;
    dA = dims(1);
    dB = dims(2);
    assert(size(rhoAB, 1) == dA * dB);
    assert(size(rhoAB, 2) == dA * dB);
    usePPT = length(ppt) > 0;
    tauPT = {};
    tauPTC = {};
    symC = [];
    reprC = [];
    if k == 1 % handle separately the PPT criterion alone, without extension
        tau = rhoAB;
        tauC = [rhoAB >= 0];
        if isequal(ppt, [])
            warning('The symmetric 1-extension without PPT constraints is trivial.');
        else
            rhoTA = reshape(rhoAB, [dB dA dB dA]);
            rhoTA = permute(rhoTA, [3 2 1 4]);
            rhoTA = reshape(rhoTA, [d d]);
            PT = {rhoTA};
            PTC = {rhoTA >= 0};
        end
    else % k > 1
        % construct the symmetric extension density matrix
        if useSym
            [~, G] = BasisSymmetricSubspace(dB, k);
            dBext = size(G, 2);
            tau = sdpvar(dA*dBext, dA*dBext, 'hermitian');
            conv = kron(eye(dA), G);
            tauFull = conv * tau * conv';
            MainCons = [tau >= 0];
        else
            dBext = dB^k;
            [~, nPi, dPi] = ProjectorSymmetricSubspace(dB, k);
            nPi = kron(eye(dA), nPi);
            tau = sdpvar(dA*dBext, dA*dBext, 'hermitian');
            tauFull = tau;
            symC = [nPi * tauFull * nPi == dPi * dPi * tauFull]; % force symmetry
        end
        % the symmetric extension is SDP
        tauC = tau >= 0;
        % the trace of the extension reproduces rhoAB
        tauAB = reshape(tauFull, [dB dB^(k-1) dA dB dB^(k-1) dA]);
        tauAB = permute(tauAB, [1 3 4 6 2 5]);
        tauAB = reshape(tauAB, [dB*dA*dB*dA dB^(k-1)*dB^(k-1)]);
        % partial trace
        tauAB = tauAB * reshape(eye(dB^(k-1)), dB^(k-1)*dB^(k-1), 1);
        tauAB = reshape(tauAB, [dB*dA dB*dA]);
        reprC = [tauAB == rhoAB];
        if usePPT
            for i = 1:length(ppt)
                k1 = ppt(i);
                k2 = k - ppt(i);
                if useSym
                    [~, G1] = BasisSymmetricSubspace(dB, k1);
                    [~, G2] = BasisSymmetricSubspace(dB, k2);
                    dB1 = size(G1, 2);
                    dB2 = size(G2, 2);
                    conv = kron(eye(dA), kron(G1, G2) \ G); % split the symmetric subspace
                    tauPT{i} = conv * tau * conv';
                    tauPT{i} = reshape(tauPT{i}, [dB1 dB2 dA dB1 dB2 dA]);
                    tauPT{i} = permute(tauPT{i}, [4 2 3 1 5 6]);
                    tauPT{i} = reshape(tauPT{i}, [dB1*dB2*dA dB1*dB2*dA]);
               
                else
                    tauPT{i} = reshape(tauFull, [dB^k1 dB^k2 dA dB^k1 dB^k2 dA]);
                    tauPT{i} = permute(tauPT{i}, [4 2 3 1 5 6]);
                    tauPT{i} = reshape(tauPT{i}, [dB^k*dA dB^k*dA]);
                end
                tauPTC{i} = [tauPT{i} >= 0];
            end
        end
    end
    Cons = [tauC
            symC
            reprC
            tauPTC{:}];
    if nargout > 1
        Info = struct;
        Info.tauC = tauC;
        Info.symC = symC;
        Info.reprC = reprC;
        Info.tauPTC = tauPTC;
        Info.tau = tau;
        Info.tauFull = tauFull;
        Info.tauPT = tauPT;
    end
end
