function GI = gouldIndex(ADJ)
    [dominant_eigval, dominant_eigvec] = epair(ADJ + diag(ones(1, width(ADJ))));
    GI = find(max(dominant_eigvec)==dominant_eigvec);
end