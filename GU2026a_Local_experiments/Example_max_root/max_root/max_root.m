% Function (6.1) from the article

function y = max_root(x)
    a = 1/4;
    y = max(sqrt(abs(x) + a) - sqrt(a));

end

