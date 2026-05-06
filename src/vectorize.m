function g = vectorize(fun, x)
    [dx, dy] = fun(x);
    g = [dx; dy];   % stack into a single column vector
end