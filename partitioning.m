clear all
n = 10;

ADJ = randi([0,1], n, n);
ADJ = ADJ - tril(ADJ,-1) + triu(ADJ,1)';
ADJ = ADJ - diag(diag(ADJ));
DEG = diag(sum(ADJ));
LAPLACIAN = DEG - ADJ;
LAP_INV = LAPLACIAN^(-1);

[V, D] = eig(LAP_INV);
eigvals = nonzeros(abs(D));
dom_eigval = max(eigvals);
i = find(dom_eigval == nonzeros(abs(D)));
eigvals(i) = 0;
dom_eigvec = V(:,i);
second_eigval = max(eigvals);
j = find(second_eigval == nonzeros(abs(D)));
Fiedler = V(:, j);

%sort_Fiedler = sort(Fiedler);

figure;
G = graph(ADJ);

subplot(2,2,1);
plot_handler = plot(G);
title("Centrality, eig");
subtitle("red node is central");
gould_index = gouldIndex(ADJ);
highlight(plot_handler, [gould_index], 'NodeColor', 'red');

subplot(2,2,2);
plot_handler = plot(G);
title("Partitioning - positive, negative, eig")
subtitle("green nodes negative; blue nodes positive");
highlight(plot_handler, [find(Fiedler < 0)], 'NodeColor', 'green');

subplot(2,2,3);
plot_handler = plot(G);
title("Partitioning - mean, eig");
subtitle("green nodes < mean; blue nodes >= mean");
split_val = mean(Fiedler);
highlight(plot_handler, [find(Fiedler < split_val)], 'NodeColor', 'green');

subplot(2,2,4);
plot_handler = plot(G);
title("Partitioning - median, eig");
subtitle("green nodes < median; blue nodes >= median");
split_val = median(sort(Fiedler));
highlight(plot_handler, [find(Fiedler < split_val)], 'NodeColor', 'green');


[dom_eigval, dom_eigvec] = epair(LAP_INV);
[second_eigval, Fiedler] = deflation(LAP_INV, dom_eigval, dom_eigvec);
figure;
G = graph(ADJ);

subplot(2,2,1);
plot_handler = plot(G);
title("Centrality, epair");
subtitle("red node is central");
gould_index = gouldIndex(ADJ);
highlight(plot_handler, [gould_index], 'NodeColor', 'red');


subplot(2,2,2);
plot_handler = plot(G);
title("Partitioning - positive, negative, epair")
subtitle("green nodes negative; blue nodes positive");
highlight(plot_handler, [find(Fiedler < 0)], 'NodeColor', 'green');

subplot(2,2,3);
plot_handler = plot(G);
title("Partitioning - mean, epair");
subtitle("green nodes < mean; blue nodes >= mean");
split_val = mean(Fiedler);
highlight(plot_handler, [find(Fiedler < split_val)], 'NodeColor', 'green');

subplot(2,2,4);
plot_handler = plot(G);
title("Partitioning - median, epair");
subtitle("green nodes < median; blue nodes >= median");
split_val = median(sort(Fiedler));
highlight(plot_handler, [find(Fiedler < split_val)], 'NodeColor', 'green');




