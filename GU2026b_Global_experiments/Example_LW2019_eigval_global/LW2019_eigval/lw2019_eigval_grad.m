function grad = lw2019_eigval_grad(x,A_tensor)

    n = size(x,1);
    [~,eigvec] = lw2019_eigval_f(x,A_tensor);

    V_mat = eigvec*eigvec';
    grad = zeros(n,1);
    for i = 1:n
        grad(i) = sum(sum(squeeze(A_tensor(i+1,:,:)) .* V_mat));
    end

end

