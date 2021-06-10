function [eigval, eigvec] = epair(A)
    N_ITER = 100;
    r = height(A);
    eigvec = ones(r,1);
    for k = 1:N_ITER
        mul = A * eigvec;
        eigvec = mul ./ norm(mul);
    end   
    eigval = (eigvec' * A * eigvec)/(eigvec' * eigvec);
end