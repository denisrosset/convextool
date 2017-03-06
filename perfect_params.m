function params = perfect_params
output_cols = [ 1  2  7 5	
                3 4 8 6	
                9 11 14 16
                10 12 13 15];
inputs_order = [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 5 5 5 5 5 5 6 6 6 6 6 6
                1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 6 6 6 6 5 5 5 5 1 2 ...
                3 4 6 5 1 2 3 4 6 5 ];
input_vectors = [ 1  0  0
                  -1  0  0
                  0  1  0
                  0 -1  0
                  0  0  1
                  0  0 -1];
params = create_params(inputs_order, output_cols, input_vectors, []);
end
