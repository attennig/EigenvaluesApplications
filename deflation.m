function [eigval, eigvec] = deflation(A, dom_eigval, dom_eigvec)
    %linar transormation that preserves the eigvals and eigvecs of M except
    %the dominant one
    B = A - dom_eigval * (dom_eigvec * dom_eigvec');
    %power method on B
    [eigval, eigvec] = epair(B);

end