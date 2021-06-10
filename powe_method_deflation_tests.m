clear all;
n_tests = 10;
DIMS = [50:50:500];
dominant_eigvec_error = zeros(1, length(DIMS));
dominant_eigval_error = zeros(1, length(DIMS));
second_eigvec_error = zeros(1, length(DIMS));
second_eigval_error = zeros(1, length(DIMS));

for n = DIMS
    for k = [1:n_tests]
        M = rand(n,n);
        [V, D] = eig(M);
        
        eigvals = nonzeros(abs(D));
        dominant_eigval = max(eigvals);
        i = find(dominant_eigval == nonzeros(abs(D)));
      
        eigvals(i) = 0;
        dominant_eigvec = V(:,i);
        second_eigval = max(eigvals);
        j = find(second_eigval == nonzeros(abs(D)));
        second_eigvec = V(:, j);
        
        % computing dominant eigenvalue/vector using power method
        [dominant_eigval_PM, dominant_eigvec_PM] = epair(M); 
        dominant_eigvec_error(n/50) =  dominant_eigvec_error(n/50) + abs(norm(dominant_eigvec) - norm(dominant_eigvec_PM))/n_tests;
        dominant_eigval_error(n/50) = dominant_eigval_error(n/50) + abs(dominant_eigval - dominant_eigval_PM)/n_tests;
        
        % computing second max eigenvalue/vector
        [second_eigval_PM, second_eigvec_PM] = deflation(M, dominant_eigval, dominant_eigvec);
        second_eigvec_error(n/50) = second_eigvec_error(n/50) + abs(norm(second_eigvec) - norm(second_eigvec_PM))/n_tests;
        second_eigval_error(n/50) = second_eigval_error(n/50) + abs(second_eigval - second_eigval_PM)/n_tests; 
    end
end

figure;
subplot(1, 2, 1);
errorbar(DIMS, dominant_eigval_error, std(dominant_eigval_error)*ones(size(DIMS)));
title("Power Method error, eigenvalues");
xlabel("Dim"); 
ylabel("Error");
subplot(1, 2, 2);
errorbar(DIMS, dominant_eigvec_error, std(dominant_eigvec_error)*ones(size(DIMS)));
title("Power Method error, eigenvectors");
xlabel("Dim"); 
ylabel("Error");
savefig("Figures/power_method_errors");

figure;
subplot(1, 2, 1);
errorbar(DIMS, second_eigval_error, std(second_eigval_error)*ones(size(DIMS)));
title("Deflation Method error, eigenvalues");
xlabel("Dim"); 
ylabel("Error");
subplot(1, 2, 2);
errorbar(DIMS, second_eigvec_error, std(second_eigvec_error)*ones(size(DIMS)));
title("Deflation Method error, eigenvectors");
xlabel("Dim"); 
ylabel("Error");
savefig("Figures/deflation_method_errors");