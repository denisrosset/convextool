function set = RandomRobustnessConeC(dims, k, ppt, useSym)
    assert(length(dims) == 2);
    dA = dims(1);
    dB = dims(2);
    d = dA*dB;
    rhoI = eye(d)/d;
    cvx_begin set sdp
    variable rho(d, d) hermitian
    variable nu nonnegative
    rho + nu * rhoI == SymmetricExtensionConeC(dims, k, ppt, useSym)
    cvx_end
    set = cvxtuple(struct('nu', nu, 'rho', rho));
end
