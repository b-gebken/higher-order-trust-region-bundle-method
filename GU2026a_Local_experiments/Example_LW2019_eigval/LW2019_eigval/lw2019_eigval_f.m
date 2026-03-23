% Function (8.7) from [Lewis, Wylie (2019), "A simple Newton method for
% local nonsmooth optimization"] 

function [eigval,eigvec] = lw2019_eigval_f(x,A_tensor)

    n = size(x,1);
    xA = squeeze(A_tensor(1,:,:));
    for i = 2:n+1
        xA = xA + x(i-1)*squeeze(A_tensor(i,:,:));
    end
    [eigvec,eigval] = eigs(xA,1,"largestreal");
    eigvec = eigvec/norm(eigvec,2);
    eigvec = eigvec/sign(eigvec(1));

end

