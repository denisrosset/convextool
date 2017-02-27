function rho = SymmetricExtensionConeC(def)
% SymmetricExtensionConeC Approximation of the cone of separable operators
%
% Formulation in the CVX standard form (= SeDuMi primal, using equalities)
% There is no support of "realification". 
%
% INPUTS
% def        Symmetric cone definition, see SymmetricExtensionDef
%
% OUTPUTS
% rho        Complex (dA*dB)x(dA*dB) matrix representing the AB system
%            The A,B basis ordering is such that a product state
%            rho = rhoA (x) rhoB = kron(rhoA, rhoB)
    
    if ~def.useSym
        error('For now, we only support the symmetric form.');
    end
    
    dims = def.dims;
    dA = dims(1);
    dB = dims(2);
    d = dA*dB;
    k = def.k;
    
    if k == 1 % handle separately the PPT criterion alone, without extension
        d = dA*dB;
        if length(ppt) == 0
            warning('The symmetric 1-extension without PPT constraint is trivial');
            cvx_begin set sdp
            variable rhoAB(d, d) semidefinite hermitian % main variable
            cvx_end
        else
            error('TODO: better implementation');
            %            cvx_begin set sdp
            %            variable rhoAB(d, d) semidefinite hermitian % main variable
            %            variable pptAB(d, d) semidefinite hermitian % partial transpose
            %            hermMatIndices = SymmetricSubspace(dA*dB, 2);
            %            hermMatIndices 
            %            rhoTA = reshape(rhoAB, [dB dA dB dA]);
            %            rhoTA = permute(rhoTA, [3 2 1 4]);
            %            rhoTA = reshape(rhoTA, [d d]);
            %            rhoTA == pptAB
            %            cvx_end
        end
        return
    end
    
    tauS = SymmetricSubspace(dB, k);
    dBsym = tauS.dim;
    [Arho AtauSym] = def.ConstraintRepresentsSym;
   
    cvx_begin set sdp
    variable rho(dB*dA, dB*dA) hermitian % main variable
    variable tauSym(dBsym*dA, dBsym*dA) hermitian % symmetric extension in symmetric basis
    rho >= 0
    tauSym >= 0
    Arho * rho(:) == AtauSym * tauSym(:)
    cvx_end
end
