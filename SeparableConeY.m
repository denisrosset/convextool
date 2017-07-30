function CONS = SeparableConeY(def, rho)
% SeparableConeC Exact/approx. formulation of the cone of separable operators
%
% INPUTS
%
% def        Separable cone exact/approx. definition, see SeparableConeDef
%
% rho        State on which the separable cone constraints apply
%            Complex (dA*dB)x(dA*dB) matrix representing the AB system
%            The A,B basis ordering is such that a product state
%            rho = rhoA (x) rhoB = kron(rhoA, rhoB)
%
% OUTPUTS
%
% CONS       Set of SDP constraints defining the separable cone formulation/approximation
    
    dims = def.dims;
    dA = dims(1);
    dB = dims(2);
    d = dA*dB;
    k = def.k;
    
    if k == 1 % handle separately the PPT criterion alone, without extension
        d = dA*dB;
        if length(def.pptCuts) == 0
            warning('The symmetric 1-extension without PPT constraint is trivial');
            CONS = SemidefiniteY(rho);           
        else
            ppt1 = sdpvar(d, d, 'hermitian', 'complex');
            [Arho Appt] = def.ConstraintPPT;
            CONS = [SemidefiniteY(rho)
                    SemidefiniteY(ppt1)
                    Arho * rho(:) == Appt * ppt1(:)];
        end
        return
    end
    
    tauS = SymmetricSubspace(dB, k);
    dBsym = tauS.dim;
    
    pptvar = cell(1, k);
    tauSym = sdpvar(dBsym*dA, dBsym*dA, 'hermitian', 'complex'); % symmetric extension in symmetric basis
    [AtauSym Arho ArhoAIdB] = def.ConstraintRepresentsSym;
    if isequal(def.approx, 'inner')
        epsN = def.innerFactor;
        CONS = [(dB*(1-epsN) * Arho + epsN * ArhoAIdB) * rho(:) == dB * AtauSym * tauSym(:)];
    else
        CONS = [Arho * rho(:) == AtauSym * tauSym(:)];
    end
    CONS = [SemidefiniteY(tauSym)
            CONS];
    for p = def.pptCuts
        k1 = p;
        k2 = k - p;
        dBsym1 = SymmetricSubspace(dB, k1).dim;
        dBsym2 = SymmetricSubspace(dB, k2).dim;
        varName = ['ppt' num2str(p)];
        ppt = sdpvar(dBsym1*dBsym2*dA, dBsym1*dBsym2*dA, 'hermitian', 'complex');
        [ApptSym AtauSym] = def.ConstraintPPTCutSym(p);
        CONS = [CONS
                SemidefiniteY(ppt)
                ApptSym * ppt(:) == AtauSym * tauSym(:)];
        pptvar{p} = ppt;
    end
end
