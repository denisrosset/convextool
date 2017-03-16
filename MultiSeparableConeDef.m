classdef MultiSeparableConeDef

    properties(SetAccess = immutable)
        
        n;            % number of subsystems
        
        dims;         % = [d1 d2 ... dn] dimensions of the subsystems composing the symmetric cone
        
        copies;       % = [k1 k2 ... kn] ki >= 1, total number of copies of each subsystem
        
        cuts;         % = m x n matrix listing the m ppt cuts, each cut is written
                      %   [p1 p2 ... pn] where pi is the number of copies of each
                      %   subsystem on which the partial transpose applies
    end

    methods
        
        function def = MultiSeparableConeDef(dims, copies, cuts)
        % MultiSeparableConeDef Defines an outer approximation of the multipartite separable cone
        %
        % INPUTS
        %
        % dims    = [d1 d2 ... dn] Dimensions of the subsystems composing the symmetric cone
        % copies  = [k1 k2 ... kn] Total number of copies of each subsystem
        % cuts    = m x n matrix listing the m ppt cuts, each cut is written
        %           [p1 p2 ... pn] where pi is the number of copies of each subsystem
        %                          on which the partial transpose applies
            assert(length(dims) == length(dims(:)));
            n = length(dims);
            dims = dims(:)';
            assert(length(copies) == length(copies(:)));
            assert(length(copies) == n);
            copies = copies(:)';
            if nargin < 3 || isequal(cuts, [])
                cuts = zeros(0, n);
            end
            assert(size(cuts, 2) == n);
            assert(all(dims >= 1));
            assert(all(copies >= 1));
            assert(all(cuts(:) >= 0));
            assert(all(all(ones(size(cuts, 1), 1)*copies - cuts >= 0)));
            def.n = n;
            def.dims = dims;
            def.copies = copies;
            def.cuts = cuts;
        end

        function [ApptSym AtauSym Dppt] = ConstraintPPTCutSym(def, p)
        % equality constraints that express that
        % tau^(partial transposes) == state
        % in their respective symmetric subspaces
        %
        % the constraint is ApptSym * pptSym(:) == AtauSym * tauSym(:)
            n = def.n;
            p = p(:)';
            q = def.copies - p;
            Stau = cell(1, n); % s is the symmetric subspace over all copies
            Sp = cell(1, n);   % p is the symmetric subspace over the ppt affected part
            Sq = cell(1, n);   % q is the symmetric subspace over the ppt nonaffected part
                               % along with their dimensions
            dtau = zeros(1, n);
            dp = zeros(1, n);
            dq = zeros(1, n);
            dpq = zeros(1, n); % dimpq(i) = dimp(i) * dimq(i) dimension of the ppt result
            for i = 1:n
                Stau{i} = SymmetricSubspace(def.dims(i), def.copies(i));
                Sp{i} = SymmetricSubspace(def.dims(i), p(i));
                Sq{i} = SymmetricSubspace(def.dims(i), q(i));
                dtau(i) = Stau{i}.dim;
                dp(i) = Sp{i}.dim;
                dq(i) = Sq{i}.dim;
                dpq(i) = dp(i) * dq(i);
            end
            % uppercase D is the total Hilbert space dimension
            Dppt = prod(dpq); % the pq state is the ppt result
            Dtau = prod(dtau); % the tau state is the symmetric state encoding the extension
            vec = @(x) x(:);
            dqp_interleave = vec([dq;dp]); % qn pn .. q2 p2 q1 p1
            Mppt = MultiIndex(dqp_interleave);
            Mtau = MultiIndex(dtau);
            MpptMat = MultiIndex([Dppt Dppt]);
            MtauMat = MultiIndex([dtau dtau]);
            nEqs = (Dppt+1)*Dppt/2;
            ApptSymT = sparse(Dppt^2, nEqs);
            AtauSymT = sparse(Dtau^2, nEqs);
            ieq = 1;
            for r = 1:Dppt
                cols = (r:Dppt)'; % we only need constraints for the upper triangle
                nEqs = length(cols);
                rows = r*ones(nEqs, 1);
                rowInd = Mppt.indToSub(rows);
                colInd = Mppt.indToSub(cols);
                rowIndTau = zeros(nEqs, n);
                colIndTau = zeros(nEqs, n);
                for i = 1:n
                    delta = 2*(i-1);
                    rowFullSub = [Sq{i}.indToSubSym(rowInd(:,delta+1)) Sp{i}.indToSubSym(colInd(:,delta+2))];
                    colFullSub = [Sq{i}.indToSubSym(colInd(:,delta+1)) Sp{i}.indToSubSym(rowInd(:,delta+2))];
                    rowIndTau(:,i) = Stau{i}.subToIndSym(rowFullSub);
                    colIndTau(:,i) = Stau{i}.subToIndSym(colFullSub);
                end                    
                blockRows = ieq:ieq+nEqs-1;
                ApptSymT(MpptMat.subToInd([rows cols]), blockRows) = eye(nEqs);
                AtauSymT(MtauMat.subToInd([rowIndTau colIndTau]), blockRows) = eye(nEqs);
                ieq = ieq + nEqs;
            end
            ApptSym = ApptSymT.';
            AtauSym = AtauSymT.';
        end

        function [AtauSym Arho Dtau] = ConstraintRepresentsSym(def)
        % equality constraints that express that
        % trace_copies(tauFull) = rho
        % with tau represented in the symmetric subspace as tauSym
        %
        % the constraint is Arho * rho(:) == AtauSym * tauSym(:)

        % comments specify the correspondence with SeparableConeDef (SCD)
            d = prod(def.dims);
            n = def.n;
            Stau = cell(1, n); % the symmetric subspace over all copies
            Mfull = cell(1, n); % the subindices for the (nonsymmetrized) copies
            Mrest = cell(1, n); % the subindices for the part we perform the partial trace over
            drest = zeros(1, n); % dimension over which we perform partial trace
            dtau = zeros(1, n); % dimension of the symmetric subspace for each subsystem
            for i = 1:n
                Stau{i} = SymmetricSubspace(def.dims(i), def.copies(i));
                Mrest{i} = MultiIndex([def.dims(i) * ones(1, def.copies(i)-1)]);
                Mfull{i} = MultiIndex([def.dims(i) * ones(1, def.copies(i))]);
                Msplit{i} = MultiIndex([def.dims(i) def.dims(i)^(def.copies(i)-1)]);
                dtau(i) = Stau{i}.dim;
                drest(i) = def.dims(i)^(def.copies(i)-1);
            end
            Mrho = MultiIndex(def.dims);
            MrhoMat = MultiIndex([d d]);
            MrestSplit = MultiIndex(drest);
            MtauMatSplit = MultiIndex([dtau dtau]);               
            Drest = prod(drest); % dimension over which we perform the partial trace
            Dtau = prod(dtau);
            traceOver = (1:Drest)';
            nEqs = (d+1)*d/2; % upper triangle dimension
            Arho = sparse(nEqs, d^2);
            AtauSymT = sparse(Dtau^2, nEqs);
            ieq = 1;
            % (r,c) move around rho(r,c)
            for r = 1:d
                rrho = Mrho.indToSub(r);
                rrest = MrestSplit.indToSub(traceOver); % indices on the "rest" part
                rsym = zeros(Drest, n);
                for i = 1:n
                    riFullSub = [rrho(i)*ones(Drest, 1) Mrest{i}.indToSub(rrest(:,i))];
                    rsym(:,i) = Stau{i}.subToIndSym(riFullSub);
                end
                for c = r:d % restrict constraints to upper triangle
                    crho = Mrho.indToSub(c);
                    crest = MrestSplit.indToSub(traceOver);
                    csym = zeros(Drest, n);
                    for i = 1:n
                        ciFullSub = [crho(i)*ones(Drest, 1) Mrest{i}.indToSub(crest(:,i))];
                        csym(:,i) = Stau{i}.subToIndSym(ciFullSub);
                    end
                    Arho(ieq, MrhoMat.subToInd([r c])) = 1;
                    tauInd = MtauMatSplit.subToInd([rsym csym]);
                    tauOnes = ones(size(tauInd));
                    col = sparse(tauInd, tauOnes, tauOnes, Dtau*Dtau, 1);
                    AtauSymT(:, ieq) = col;
                    ieq = ieq + 1;
                end
            end
            AtauSym = AtauSymT.';
        end
    end

end
