function rhoAB = SymmetricExtensionConeC(dims, k, ppt, useSym)
% SymmetricExtensionConeC Approximation of the cone of separable operators
%
% Formulation in the CVX standard form (= SeDuMi primal, using equalities)
% To gain advantage of this formulation, set sdpsettings('dualize', 1)
% There is no support of "realification". This formulation is quite
% inefficient for `useSym = 0`.
%
% INPUTS
% dims       = [dA dB]  Dimension of subsystems A and B
% k          Number of copies of subsystem B in the extension
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
%
% OUTPUTS
% rhoAB      Complex (dA*dB)x(dA*dB) matrix representing the AB system
%            The A,B basis ordering is such that a product state
%            rhoAB = rhoA (x) rhoB = kron(rhoA, rhoB)

    % process parameters
    assert(length(dims) == 2);
    dA = dims(1);
    dB = dims(2);
    if nargin < 3
        ppt = [];
    end
    if nargin < 4 || isequal(useSym, [])
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

    % start the convex cone definition
    cvx_begin set sdp
    variable rhoAB(dA*dB, dA*dB) hermitian % main variable

    if useSym
        [~, G] = BasisSymmetricSubspace(dB, k);
        dBext = size(G, 2);
        variable tau(dA*dBext, dA*dBext) hermitian % variable
        conv = kron(eye(dA), G);
        tauFull = conv * tau * conv';
        tau >= 0; % semidefinite positive
    else
        dBext = dB^k;
        [~, nPi, dPi] = ProjectorSymmetricSubspace(dB, k);
        nPi = kron(eye(dA), nPi);
        variable tauFull(dA*dBext, dA*dBext) hermitian % variable       
        nPi * tauFull * nPi == dPi * dPi * tauFull; % symmetry
        tauFull >= 0; % semidefinite positive
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
                tauPPTvar >= 0;
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
                tauPPTvar >= 0;
            end
        end
    end
    cvx_end
end
