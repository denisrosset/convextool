function [consistentP_a_b_i, lsq_error] = make_prob_consistent(P_a_b_i, params)
    comp_a_b_i = sdpvar(params.nA, params.nB, params.nI);
    CONS = [];
    Piab = cell(params.nA, params.nB);
    sumPiR = zeros(8, 8);
    for a = 1:params.nA
        for b = 1:params.nB
            [Pi, PiR] = hermitian2(sdpvar(16, 1));
            Piab{a, b} = Pi;
            sumPiR = sumPiR + PiR;
            CONS = [CONS
                    PiR >= 0];
            for inp = 1:params.nI
                comp_a_b_i(a,b,inp) = real(trace(Piab{a, b} * params.input_states{inp}));
            end
        end
    end
    CONS = [CONS
            eye(8) - sumPiR >= 0];
    obj = norm(comp_a_b_i(:) - P_a_b_i(:), 2);
    solvesdp(CONS, obj, sdpsettings('verbose', 1));
    consistentP_a_b_i = double(comp_a_b_i);
    lsq_error = double(obj);
end
