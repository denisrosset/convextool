function [rho CONS] = MultiSeparableConeY(def)
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
    [CR_AtauSym CR_Arho Dtau] = def.ConstraintRepresentsSym(false);
    nPPT = size(def.cuts, 1);
    pptvar = cell(1, nPPT);
    tauSym = sdpvar(Dtau, Dtau, 'hermitian', 'complex'); % symmetric extension in symmetric basis
    drho = prod(def.dims);
    rho = reshape(CR_AtauSym * tauSym(:), [drho drho]).'; % because CR_Arho is eye(d) transposed
       
    CONS = [SemidefiniteY(tauSym)];   
        
    for p = 1:nPPT
        [ApptSym AtauSym Dppt] = def.ConstraintPPTCutSym(def.cuts(p,:), false);
        ppt = reshape(AtauSym * tauSym(:), [Dppt Dppt]).'; % because ApptSym is eye(d) transposed
        CONS = [CONS
                SemidefiniteY(ppt)];
        pptvar{p} = ppt;
    end
end
