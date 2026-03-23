function grad = grad_halfhalf(x)

    A = diag([1,0,1,0,1,0,1,0]);
    B = diag(1./(1:8).^2);

    grad = (x'*A*x)^(-1/2)*A*x + 2*B*x;

end

