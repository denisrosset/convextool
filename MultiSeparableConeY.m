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
    [CR_AtauSym CR_Arho Dtau] = def.ConstraintRepresentsSym;
    nPPT = size(def.cuts, 1);
    pptvar = cell(1, nPPT);
    
    tauSym = sdpvar(Dtau, Dtau, 'hermitian', 'complex'); % symmetric extension in symmetric basis
    drho = prod(def.dims);
    rho = sdpvar(drho, drho, 'hermitian', 'complex');
    for cstr = 1:size(CR_Arho, 1)
        ind = find(CR_Arho(cstr, :));
        assert(length(ind) == 1);
        [r c] = ind2sub([drho drho], ind);
        rho(r,c) = CR_AtauSym(cstr,:)*tauSym(:);
        if r ~= c
            rho(c,r) = conj(CR_AtauSym(cstr,:)*tauSym(:));
        end
    end
   
    CONS = [SemidefiniteY(tauSym)];   
        
    for p = 1:nPPT
        [ApptSym AtauSym Dppt] = def.ConstraintPPTCutSym(def.cuts(p,:));
        ppt = sdpvar(Dppt, Dppt, 'hermitian', 'complex');
        for cstr = 1:size(ApptSym, 1)
            ind = find(ApptSym(cstr, :));
            assert(length(ind) == 1);
            [r c] = ind2sub([Dppt Dppt], ind);
            ppt(r,c) = AtauSym(cstr,:)*tauSym(:);
            if r ~= c
                ppt(c,r) = conj(AtauSym(cstr,:)*tauSym(:));
            end
        end
        CONS = [CONS
                SemidefiniteY(ppt)];
        pptvar{p} = ppt;
    end
end
