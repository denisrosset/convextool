function [Cons MainCons PPTCons] = SymmetricExtensionConePrimal(rhoAB, dA, dB, k, ppt, useSym)
% SymmetricExtensionConePrimal Compute SDP constraints corresponding to symmetric extension cones
%
% Formulation in the YALMIP primal canonical form (= using equalities)
% To gain advantage of this method, set sdpsettings('dualize', 1)
% 
%
% INPUTS
% rhoAB      (dA*dB)x(dA*dB) matrix representing the AB system
%            The A,B basis ordering is such that a product state
%            rhoAB = rhoA (x) rhoB = kron(rhoA, rhoB)
% dA         Dimension of subsystem A
% dB         Dimension of subsystem B
% k          Number of copies of subsystem B
% ppt        PPT constraints to add. Can be of the following form:
%            = []            no PPT constraints
%            = [t_1 ... t_n] with 1 <= t_j <= k (duplicates ignored)
%                            for each t_j, adds a PPT constraint such that a number t_j of copies
%                            of B is transposed
%            = 'doherty'     equivalent to [1 2 ... k], PPT conditions in the original 2004 Doherty paper
%            = 'navascues'   equivalent to [ceil(k/2)], PPT condition in the 2009 Navascues paper,
%                                                       also the way it is implemented in QETLAB
%            Default: [] (do not use PPT constraints)
% useSym     Whether to use the symmetric subspace (default: 1)
    assert(size(rhoAB, 1) == dA * dB);
    assert(size(rhoAB, 2) == dA * dB);
    if nargin < 5
        ppt = [];
    end
    if nargin < 6 || isequal(useSym, [])
        useSym = true;
    end
    if isequal(ppt, 'doherty')
        ppt = 1:k;
    elseif isequal(ppt, 'navascues')
        ppt = ceil(k/2);
    else
        assert(all(ppt >= 1));
        assert(all(ppt <= k));
    end
    usePPT = length(ppt) > 0;
    if usePPT
        ppt = unique(ppt);
    end
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
        tauFull = sdpvar(dA*dBext, dA*dBext, 'hermitian');
        MainCons = [nPi * tauFull * nPi == dPi * dPi * tauFull % force symmetry
                    tauFull >= 0];
    end
    tauAB = reshape(tauFull, [dB dB^(k-1) dA*fieldDim dB dB^(k-1) dA*fieldDim]);
    tauAB = permute(tauAB, [1 3 4 6 2 5]);
    tauAB = reshape(tauAB, [dB*dA*fieldDim*dB*dA*fieldDim dB^(k-1)*dB^(k-1)]);
    % partial trace
    tauAB = tauAB * reshape(eye(dB^(k-1)), dB^(k-1)*dB^(k-1), 1);
    tauAB = reshape(tauAB, [dB*dA*fieldDim dB*dA*fieldDim]);
    if realify
        rhoRealPart = real(rhoAB);
        rhoImagPart = imag(rhoAB);
        [tauRealPart tauImagPart] = RealifiedParts(tauAB);
        MainCons = [MainCons
                    rhoRealPart == tauRealPart
                    rhoImagPart == tauImagPart];
        %        MainCons = [MainCons
        %                    tauAB == ComplexToReal(rhoAB)];
    else
        MainCons = [MainCons
                    tauAB == rhoAB];
    end
    PPTCons = [];
    if usePPT
        if useSym
            for i = 1:length(ppt)
                k1 = ppt(i);
                k2 = k - ppt(i);
                [~, G1] = BasisSymmetricSubspace(dB, k1);
                [~, G2] = BasisSymmetricSubspace(dB, k2);
                dB1 = size(G1, 2);
                dB2 = size(G2, 2);
                conv = kron(kron(eye(dA), kron(G1, G2) \ G), fieldId); % split the symmetric subspace
                tauPPT = conv * tau * conv';
                tauPPT = reshape(tauPPT, [dB1 dB2 dA*fieldDim dB1 dB2 dA*fieldDim]);
                tauPPT = permute(tauPPT, [4 2 3 1 5 6]);
                tauPPT = reshape(tauPPT, [dB1*dB2*dA*fieldDim dB1*dB2*dA*fieldDim]);
                PPTCons = [PPTCons
                           tauPPT >= 0];
            end
        else
            for i = 1:length(ppt)
                k1 = ppt(i);
                k2 = k - ppt(i);
                tauPPT = reshape(tauFull, [dB^k1 dB^k2 dA*fieldDim dB^k1 dB^k2 dA*fieldDim]);
                tauPPT = permute(tauPPT, [4 2 3 1 5 6]);
                tauPPT = reshape(tauPPT, [dB^k*dA*fieldDim dB^k*dA*fieldDim]);
                PPTCons = [PPTCons
                           tauPPT >= 0];
            end
        end
    end
    Cons = [MainCons
            PPTCons];
end