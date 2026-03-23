% Function (8.4) from [Lewis, Wylie (2019), "A simple Newton method for
% local nonsmooth optimization"] 

function [val,I] = lw2019_84_f(x,g_arr,H_cell,c_arr)

    k = size(g_arr,2);

    tmp = zeros(k,1);
    for i = 1:k
        tmp(i) = g_arr(:,i)'*x + 1/2 * x'*H_cell{i}*x + c_arr(i)/24 * norm(x,2)^4;
    end
    
    [val,I] = max(tmp);
end

