)n = 10;

ADJ = randi([0,1], n, n);
ADJ = ADJ - tril(ADJ,-1) + triu(ADJ,1)';
ADJ = ADJ - diag(diag(ADJ));
gould_index = gouldIndex(ADJ);
G = graph(ADJ);
plot_handler = plot(G);

highlight(plot_handler, [gould_index], 'NodeColor', 'red');