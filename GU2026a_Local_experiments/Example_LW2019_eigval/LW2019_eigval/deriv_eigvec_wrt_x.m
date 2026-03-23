function jac = deriv_eigvec_wrt_x(xV,xD,A_tensor)

n = size(A_tensor,1)-1;
m = size(xV,1);

jac = zeros(m,n);
for k = 1:m
    for l = 1:m
        dv = deriv_eigvec_wrt_matr_comp(xV,xD,k,l);
        jac = jac + dv * A_tensor(2:n+1,k,l)';
    end
end

end

