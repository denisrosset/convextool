function [Cons Info] = SymmetricExtensionConeDualY(coeffs, def)
% SymmetricExtensionConeY Compute SDP constraints corresponding to symmetric extension cones
%
% Formulation in the YALMIP dual canonical form (= using inequalities); to recover the dual
% variables associated with constraints, set 'toReal' = 1; otherwise, YALMIP performs the
% real transform itself, and cannot recover them.
%
% Note that the input density matrix is specified using a real decomposition over an Hermitian
% basis, avoiding the use of complex numbers when those are undesirable.
%
% INPUTS
% coeffs     (dA^2)x(dB^2) matrix of type double/sdpvar containing the coefficients
%            of the Hermitian operator to constrain in the symmetric extension cone
%            H = sum_ij coeffs(i,j) * kron(FA{i}, FB{j})
%            where FA, FB are given by GeneralizedGellMann
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
% Info.tau =   Density matrix of the symmetric extension
%
%              For 'useSym' = 1
% Info.tau     Density matrix of the symmetric extension (reduced sym. basis)
%
%              For length(def.ppt) > 1
% Info.tauPT   Cell array of matrices of partial transposes
% Info.tauPTC  Cell array of PPT constraints PT{i} >= 0
%
%              For k > 1
% Info.reprC   The constraint that trace_{B2..Bk}(tau) == rhoAB
    dims = def.dims;
    k = def.k;
    ppt = def.ppt;
    useSym = def.useSym;
    toReal = def.toReal;
    if toReal
        fieldDim = 2; % 1 if complex, 2 if real
    else
        fieldDim = 1;
    end
    dA = dims(1);
    dB = dims(2);
    assert(size(coeffs, 1) == dA * dA);
    assert(size(coeffs, 2) == dB * dB);
    [FA DA indPTA] = GeneralizedGellMann(dA);
    [FB DB indPTB] = GeneralizedGellMann(dB);
    % indices of the rows/columns to preserve, because the
    % matrix we consider has support and range in the symmetric
    % subspace
    %
    % we consider the same for the PPT matrices
    symIndices = [];
    symIndicesPPT = cell(1, k);
    tauPT = {};
    tauPTC = {};
    usePPT = length(ppt) > 0;
    if useSym
        indicesB = SymmetricCanonicalIndices(dB, k);
        symIndices = indicesB(:)';
        dBext = length(symIndices);
        dBPPT = zeros(1, k);
        if usePPT
            % order is B1 ... Bk1 Bk1+1 ... Bk1+k2
            % but warning: kron reverses the order compared to ind2sub
            for k1 = 1:k
                k2 = k - k1;
                % thus indices1 has factor
                indices1 = (SymmetricCanonicalIndices(dB, k1) - 1) * dB^k2;
                indices2 = SymmetricCanonicalIndices(dB, k2);
                indicesB = bsxfun(@plus, indices1', indices2);
                symIndicesPPT{k1} = indicesB(:)';
                dBPPT(k1) = length(symIndicesPPT{k1});
            end
        end
    else
        dBext = dB^k;
        dBPPT = ones(1, k) * dB^k;
    end
    % order is A B1 ... Bk
    % The formulation below has been optimized to reduce preprocessing time
    % Look at commit e8ed3890 for the more explicit, slower, version
    sdpdim = dA^2*dBext^2;
    if usePPT
        sdpdim = sdpdim + sum(dA^2*dBPPT.^2);
    end
    fixed = zeros(sdpdim, 1);
    basis = sparse(sdpdim, 0);
    binds = SymmetricCanonicalIndices(dB^2, k);
    for bind = binds
        bs = cell(1, k);
        [bs{:}] = ind2sub(dB^2*ones(1, k), bind);
        bs = cell2mat(bs);
        switch sum(bs > 1)
          case 0 % Alice marginal
            b = 1;
          case 1
            b = bs(find(bs > 1));
          otherwise
            b = 0;
        end
        allbs = unique(perms(bs), 'rows');
        B = sparse(dB^k, dB^k);
        for r = 1:size(allbs, 1)
            El = 1;
            for c = 1:k
                El = kron(El, FB{allbs(r, c)});
            end
            B = B + El;
        end
        for a = 1:dA^2
            A = FA{a};
            colVec = @(x) x(:);
            if useSym
                current = colVec(kron(A, B(symIndices, symIndices)));
            else
                current = colVec(kron(A, B));
            end
            if usePPT
                PT = full(B); % partial transpositions are added one by one
                for j = 1:k
                    d1 = dB^(j-1);
                    d2 = dB;
                    d3 = dB^(k-j);
                    PTsplit = reshape(PT, [d1 d2 d3 d1 d2 d3]);
                    PTsplit = permute(PTsplit, [1 5 3 4 2 6]);
                    PT = reshape(PTsplit, [dB^k dB^k]);
                    if useSym
                        current = [current; colVec(kron(A, PT(symIndicesPPT{j}, symIndicesPPT{j})))];
                    else
                        current = [current; colVec(kron(A, PT))];
                    end
                end
            end
            if b == 0
                basis(:, end + 1) = current;
            else
                fixed = fixed + coeffs(a, b) * current;
            end
        end
    end
    %    basis = double(orth(sym(full(basis)), 'skipnormalization'));
    allinall = fixed + basis * sdpvar(size(basis, 2), 1);
    dim = dA^2*dBext^2;
    tau = reshape(allinall(1:dim), dA*dBext, dA*dBext);
    tau = (tau + tau')/2; % force Hermitian (is it really needed?)
    allinall = allinall(dim+1:end);
    if usePPT
        for j = 1:k
            dim = dA^2*dBPPT(j)^2;
            PPT{j} = reshape(allinall(1:dim), dA*dBPPT(j), dA*dBPPT(j));
            PPT{j} = (PPT{j} + PPT{j}')/2;
            allinall = allinall(dim+1:end);
        end
    end
    if toReal
        tau = ComplexToReal(tau);
    end
    tauC = [tau >= 0];
    if usePPT
        for i = 1:length(ppt)
            if toReal
                tauPT{i} = ComplexToReal(PPT{ppt(i)});
            else
                tauPT{i} = PPT{ppt(i)};
            end
            tauPTC{i} = [tauPT{i} >= 0];
        end
    end
    Cons = [tauC
            tauPTC{:}];
    if nargout > 1
        Info = struct;
        Info.tau = tau;
        Info.tauPT = tauPT;
        Info.tauC = tauC;
        Info.tauPTC = tauPTC;
    end
end
