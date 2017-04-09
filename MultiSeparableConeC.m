function rho = MultiSeparableConeC(def)
% MultiSeparableConeC Doherty approximation of the cone of multipartite separable operators
%
% Source: Doherty et al., DOI 10.1103/PhysRevA.71.032333
%
% INPUTS
% def        Multipartite separable cone approx. definition, see MultiSeparableConeDef
%
% OUTPUTS
% rho        Complex (dA*dB*...)x(dA*dB*...) matrix representing the AB... system
%            The A,B,... basis ordering is such that a product state
%            rho = rhoA (x) rhoB (x) ... = kron(..., kron(rhoB, rhoA))
%            So CURRENTLY THE OPPOSITE OF THE SEPARABLE CONE IMPLEMENTATION
    d = prod(def.dims);    
    [CR_AtauSym CR_Arho Dtau] = def.ConstraintRepresentsSym;
    nPPT = size(def.cuts, 1);
    pptvar = cell(1, nPPT);
    cvx_begin set sdp
        variable rho(d, d) hermitian % main variable
        variable tauSym(Dtau, Dtau) hermitian  % symmetric extension in symmetric basis

        rho >= 0
        tauSym >= 0
        CR_Arho * rho(:) == CR_AtauSym * tauSym(:)
        
        for p = 1:nPPT
            [ApptSym AtauSym Dppt] = def.ConstraintPPTCutSym(def.cuts(p,:));
            varName = ['ppt' num2str(p)];
            varDecl = [varName '(Dppt, Dppt)'];
            variable(varDecl, 'hermitian');
            pptvar{p} = eval(varName);
            pptvar{p} >= 0
            ppt = pptvar{p};
            ApptSym * ppt(:) == AtauSym * tauSym(:)
        end
    cvx_end
end
