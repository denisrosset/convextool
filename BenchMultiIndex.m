n = 10;
dims = randi(8, 1, n);
m = 10000;
ind = uint64(randi(prod(dims), m, 1));
mi_unopt = MultiIndex(dims, 0);
mi_opt = MultiIndex(dims, 1);
tic
for i = 1:20
    mi_unopt.indToSub(ind);
end
toc
tic
for i = 1:30
    mi_opt.indToSub(ind);
end
toc
