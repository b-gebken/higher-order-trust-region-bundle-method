function hess = hess_max_root(x)
    a = 1/4;
    [~,I] = max(sqrt(abs(x) + a) - sqrt(a));

    hess = zeros(size(x,1));

    hess(I,I) = -1/4*(abs(x(I)) + a)^(-3/2);
end

