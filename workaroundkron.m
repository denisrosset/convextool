function res = workaroundkron(x, y)
    x1 = size(x, 1);
    x2 = size(x, 2);
    y1 = size(y, 1);
    y2 = size(y, 2);
    res = y(:) * x(:)';
    res = reshape(res, [y1 y2 x1 x2]);
    res = permute(res, [1 3 2 4]);
    res = reshape(res, [y1*x1 y2*x2]);
end
