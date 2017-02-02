function [Cons MainCons PPTCons] = SymmetricExtensionCone(coeffs, k, useSym, usePPT)
% SymmetricExtensionCone Compute SDP constraints corresponding to symmetric extension cones
%
% INPUTS
% coeffs     (dA^2)x(dB^2) matrix of type double/sdpvar containing the coefficients
%            of the Hermitian operator to constrain in the symmetric extension cone
%            H = sum_ij coeffs(i,j) * kron(FA{i}, FB{j})
%            where FA, FB are given by GeneralizedGellMann
% k          Number of copies
% useSym     Whether to use the symmetric subspace
% usePPT     Whether to use the PPT constraints
    dA = sqrt(size(coeffs, 1));
    dB = sqrt(size(coeffs, 2));
    [FA DA indPTA] = GeneralizedGellMann(dA);
    [FB DB indPTB] = GeneralizedGellMann(dB);
    factor = DB(1)^(k-1); % not used (because we define a cone), but the symmetric extension has
                          % larger trace by this factor (because the basis is not normalized)
    [~, B01] = BasisSymmetricSubspace(dB, k);
    DupB = B01'*B01;
    
    
    % indices of the rows/columns to preserve, because the
    % matrix we consider has support and range in the symmetric
    % subspace
    %
    % we consider the same for the PPT matrices
    symIndices = [];
    symIndicesPPT = cell(1, k);
    if useSym
        symIndices = bsxfun(@plus, (0:dA-1)'*dB^k, SymmetricCanonicalIndices(dB, k));
        symIndices = symIndices(:)';
        S = sparse(length(symIndices)*2, length(symIndices)*2);
        if usePPT
            for k1 = 1:k
                k2 = k - k1;
                for ind = 1:dA*dB^k
                    bs1 = cell(1, k1);
                    bs2 = cell(1, k2);
                    [bs1{:} bs2{:} a] = ind2sub([dB*ones(1, k) dA], ind); % kron inverts indices
                    bs1 = cell2mat(bs1);
                    bs2 = cell2mat(bs2);
                    if all(bs1(2:end) - bs1(1:end-1) >= 0) && all(bs2(2:end) - bs2(1:end-1) >= 0)
                        symIndicesPPT{k1} = [symIndicesPPT{k1} ind];
                    end
                end
                PPT{k1} = sparse(length(symIndicesPPT{k1})*2, length(symIndicesPPT{k1})*2);
            end
        end
    else
        d = dA*dB^k;
        symIndices = 1:d;
        S = sparse(d*2, d*2);
        PPT = cell(1, k);
        for k1 = 1:k
            PPT{k1} = sparse(d*2, d*2);
        end
    end
    coeffsExt = cell(dA^2*(dB^2)^k, 1);
    range = sparse(length(symIndices), 0);
    for a = 1:dA^2
        A = FA{a};
        for bind = 1:(dB^2)^k
            bs = cell(1, k);
            [bs{:}] = ind2sub(dB^2*ones(1, k), bind);
            bs = cell2mat(bs);
            if all(bs(2:end) - bs(1:end-1) >= 0) % is increasing
                switch sum(bs > 1)
                  case 0 % Alice marginal
                    coeff = coeffs(a, 1);
                  case 1
                    coeff = coeffs(a, bs(find(bs > 1)));
                  otherwise
                    coeff = sdpvar;
                end
                allbs = unique(perms(bs), 'rows');
                for r = 1:size(allbs, 1)
                    B = 1;
                    signPPT = ones(1, k);
                    for c = 1:k
                        if indPTB(allbs(r, c))
                            signPPT(1:c) = -signPPT(1:c);
                        end
                        B = kron(B, FB{allbs(r, c)});
                    end
                    M = kron(A, B);
                    if useSym
                        S = S + coeff * ComplexToReal(M(symIndices, symIndices));
                        if usePPT
                            for k1 = 1:k
                                PPT{k1} = PPT{k1} + signPPT(k1) * coeff * ...
                                          ComplexToReal(M(symIndicesPPT{k1}, symIndicesPPT{k1}));
                            end
                        end
                    else
                        S = S + coeff * ComplexToReal(M);
                        if usePPT
                            for k1 = 1:k
                                PPT{k1} = PPT{k1} + signPPT(k1) * coeff * ComplexToReal(M);
                            end
                        end
                    end
                end
            end
        end
    end
    MainCons = [S >= 0];
    % PPT constraints
    PPTCons = [];
    if usePPT
        for k1 = 1:k
            PPTCons = [PPTCons; PPT{k1} >= 0];
        end
    end
    Cons = [MainCons
            PPTCons];
end
