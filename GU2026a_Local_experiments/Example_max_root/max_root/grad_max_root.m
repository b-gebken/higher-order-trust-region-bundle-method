function grad = grad_max_root(x)
    a = 1/4;
    [~,I] = max(sqrt(abs(x) + a) - sqrt(a));

    grad = zeros(size(x,1),1);

    s = sign(x(I));
    if(s == 0)
        s = 1;
    end

    grad(I) = 1/2*(abs(x(I)) + a)^(-1/2)*s;
end

