function [eigval, eigvec] = epair(A)
    N_ITER = 20000;
    r = height(A);
    %eigvec = ones(r,1);
    %for k = 1:N_ITER
    %    mul = A * eigvec;
    %    eigvec = mul ./ norm(mul);
    %end   
    %eigval = (eigvec' * A * eigvec)/(eigvec' * eigvec);

    x_k = ones(r,1);
    for k=1:N_ITER
        y_k = A * x_k;
        x_k = y_k ./ norm(y_k, Inf);
    end
    eigvec = x_k;
    eigval = norm(y_k, Inf);
end