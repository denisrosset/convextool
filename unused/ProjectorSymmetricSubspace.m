function [Pi nPi dPi] = ProjectorSymmetricSubspace(d, k)
    [B B01] = BasisSymmetricSubspace(d, k);
    S = sum(B01, 1);
    dPi = lcm_rep(full(S));
    nPi = sparse(B01 * diag(dPi./S) * B01');
    Pi = nPi ./ dPi;
    function res = lcm_rep(x)
        res = 1;
        for i = 1:length(x)
            res = lcm(x(i), res);
        end
    end
end