function CONS = MultiSeparableConeY(def, rho)
% MultiSeparableConeC Doherty approximation of the cone of multipartite separable operators
%
% Source: Doherty et al., DOI 10.1103/PhysRevA.71.032333
%
% INPUTS
%
% def        Multipartite separable cone approx. definition, see MultiSeparableConeDef
%
% rho        State on which the separable cone constraints apply
%            Complex (dA*dB*...)x(dA*dB*...) matrix representing the AB... system
%            The A,B,... basis ordering is such that a product state
%            rho = rhoA (x) rhoB (x) ... = kron(rhoA, kron(rhoB, ...))
%
% OUTPUTS
%
% CONS       Set of SDP constraints defining the separable cone formulation/approximation
       
    d = prod(def.dims);    
    [CR_AtauSym CR_Arho Dtau] = def.ConstraintRepresentsSym;
    nPPT = size(def.cuts, 1);
    pptvar = cell(1, nPPT);
    
    tauSym = sdpvar(Dtau, Dtau, 'hermitian', 'complex'); % symmetric extension in symmetric basis
    CONS = [SemidefiniteY(rho)
            SemidefiniteY(tauSym)
            CR_Arho * rho(:) == CR_AtauSym * tauSym(:)];
        
    for p = 1:nPPT
            [ApptSym AtauSym Dppt] = def.ConstraintPPTCutSym(def.cuts(p,:));
            ppt = sdpvar(Dppt, Dppt, 'hermitian', 'complex');
            CONS = [CONS
                    SemidefiniteY(ppt)
                    ApptSym * ppt(:) == AtauSym * tauSym(:)];
            pptvar{p} = ppt;
    end
end
