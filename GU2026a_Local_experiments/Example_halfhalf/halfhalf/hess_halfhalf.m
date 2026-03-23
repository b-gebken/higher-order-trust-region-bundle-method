function hess = hess_halfhalf(x)

    A = diag([1,0,1,0,1,0,1,0]);
    B = diag(1./(1:8).^2);

    hess = (-(x'*A*x)^(-3/2)*A*x)*(A*x)' + (x'*A*x)^(-1/2)*A + 2*B;

end

