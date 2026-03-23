function deriv = deriv_eigvec_wrt_matr_comp(V,D,k,l)

n = size(V,1);
deriv = zeros(n,1);
for j_se = 1:n
    if(j_se ~= 1)
        deriv = deriv + (V(k,j_se) * V(l,1))/(D(1,1) - D(j_se,j_se)) * V(:,j_se);
    end
end

end

