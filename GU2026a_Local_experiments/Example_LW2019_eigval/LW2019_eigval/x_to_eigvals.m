function eigvals = x_to_eigvals(x,A_tensor)

    n = size(x,1);
    xA = squeeze(A_tensor(1,:,:));
    for i = 2:n+1
        xA = xA + x(i-1)*squeeze(A_tensor(i,:,:));
    end
    D = eig(xA);
    [~,I] = sort(-D);
    eigvals = D(I);

end

