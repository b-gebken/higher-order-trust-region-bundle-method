function hess = lw2019_eigval_hess(x,A_tensor)

n = size(x,1);
m = size(A_tensor,2);

xA = squeeze(A_tensor(1,:,:));
for i = 2:n+1
    xA = xA + x(i-1) * squeeze(A_tensor(i,:,:));
end

[xV,xD] = eig(xA);
[~,ind] = sort(diag(-xD));
xD = xD(ind,ind);
xV = xV(:,ind);
xV = sign(xV(1,:)) .* xV; % Normalize first entry to positive

jac = deriv_eigvec_wrt_x(xV,xD,A_tensor);
[~,v] = lw2019_eigval_f(x,A_tensor);

hess = zeros(n);
for k = 1:m
    for l = 1:m
        hess = hess + A_tensor(2:n+1,k,l) .* (jac(k,:) * v(l) + v(k) * jac(l,:));
    end
end



end

