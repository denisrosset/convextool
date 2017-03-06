function params = create_params(input_order, output_order, input_vectors, input_mask)
    nA = size(output_order, 1);
    nB = size(output_order, 2);
    nI = size(input_order, 2);
    nS = size(input_vectors, 1);
    if ~isequal(input_mask, [])
        input_order = input_order(:, input_mask ~= 0); 
    end
    assert(size(input_vectors, 2) == 3);
    assert(size(input_order, 1) == 2);
    assert(min(input_order(:)) >= 1);
    assert(max(input_order(:)) <= nS);
    input_states = cell(nI, 1);
    for i = 1:nI
        x = input_order(1, i);
        y = input_order(2, i);
        vec_x = input_vectors(x, :);
        vec_y = input_vectors(y, :);
        tau_x = qubit_state(vec_x);
        tau_y = qubit_state(vec_y);
        input_states{i} = kron(tau_x, tau_y);
    end
    params = struct('nA', nA, 'nB', nB, 'nI', nI, 'nS', nS);
    params.input_states = input_states;
    params.output_order = output_order;
end
