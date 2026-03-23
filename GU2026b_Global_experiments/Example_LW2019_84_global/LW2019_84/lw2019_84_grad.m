function grad = lw2019_84_grad(x,g_arr,H_cell,c_arr)
    
    [~,I] = lw2019_84_f(x,g_arr,H_cell,c_arr);

    grad = g_arr(:,I) + H_cell{I}*x + c_arr(I)/24 * 4*norm(x,2)^2 * x;

end

