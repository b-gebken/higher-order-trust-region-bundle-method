% Evaluates the Taylor expansion corresponding to the derivatives in
% deriv_arr with center y at z

function vals = eval_Taylor(z,deriv_arr,y)

    q = size(deriv_arr,1) - 1;

    vals = deriv_arr' * ((1./(factorial(0:q)))' .* (((z - y)').^(0:q))');
end

