clear all

n = 15;

ADJ = randi([0,1], n, n);
ADJ = ADJ - tril(ADJ,-1) + triu(ADJ,1)';
ADJ = ADJ - diag(diag(ADJ));
DEG = diag(sum(ADJ));
LAPLACIAN = DEG - ADJ;
LAP_INV = LAPLACIAN^(-1);
[dom_eigval, dom_eigvec] = epair(LAP_INV);
[second_eigval, Fiedler] = deflation(LAP_INV, dom_eigval, dom_eigvec);

%sort_Fiedler = sort(Fiedler);

figure;
G = graph(ADJ);


subplot(2,2,1);
plot_handler = plot(G);
title("Centrality");
subtitle("red node is central");
gould_index = gouldIndex(ADJ);
highlight(plot_handler, [gould_index], 'NodeColor', 'red');


subplot(2,2,2);
plot_handler = plot(G);
title("Partitioning - positive, negative")
subtitle("green nodes negative; blue nodes positive");
highlight(plot_handler, [find(Fiedler < 0)], 'NodeColor', 'green');

subplot(2,2,3);
plot_handler = plot(G);
title("Partitioning - mean");
subtitle("green nodes < mean; blue nodes >= mean");
split_val = mean(Fiedler);
highlight(plot_handler, [find(Fiedler < split_val)], 'NodeColor', 'green');

subplot(2,2,4);
plot_handler = plot(G);
title("Partitioning - median");
subtitle("green nodes < median; blue nodes >= median");
split_val = median(sort(Fiedler));
highlight(plot_handler, [find(Fiedler < split_val)], 'NodeColor', 'green');


