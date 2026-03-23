function hess = lw2019_84_hess(x,g_arr,H_cell,c_arr)
    
    n = size(x,1);

    [~,I] = lw2019_84_f(x,g_arr,H_cell,c_arr);

    hess = H_cell{I} + c_arr(I)/24 * (8*(x*x') + diag(4*norm(x,2)^2 * ones(1,n)));

end

