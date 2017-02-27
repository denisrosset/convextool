for s1 = 1:4
    for s2 = 1:4
        for t1 = 1:4
            for t2 = 1:4
                x = rand(s1, s2); 
                y = rand(t1, t2);
                assert(all(all(workaround_kron(x, y) - kron(x, y) == 0)));
            end
        end
    end
end
