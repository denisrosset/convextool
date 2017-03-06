function P_a_b_i = estimate_prob(clicks_i_ab, params)
% ESTIMATE_PROB - estimate probablities from event counts
    assert(size(clicks_i_ab, 1) == params.nI);
    clicks_i_ab = clicks_i_ab(:, params.output_order(:));
    clicks_i_ab_abi = permute(reshape(clicks_i_ab, params.nI, params.nA, params.nB), [2 3 1]);
    sum_ab = sum(sum(clicks_i_ab_abi, 2), 1);
    sum_max = max(sum_ab);
    P_a_b_i = clicks_i_ab_abi / sum_max;
end
