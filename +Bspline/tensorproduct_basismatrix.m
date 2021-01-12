function [T,knots1,knots2] = tensorproduct_basismatrix(X, nr_splines, ll, knot_type)
%% 
% Generate the tensor-product B-spline basis given the data X!
% 
% Parameters
% ----------
% X : matrix            - Input data of shape (n_samples, 2).
% nr_splines : array    - Number of parameters for dimension 1 and 2.
% ll : array            - Order of the B-spline basis functions for 
%                         each dimension.
% knot_type : array     - Knot type of the B-spline basis functions 
%                         for each dimension 
%
% Returns:
% --------
% T  : matrix       - Tensor-product B-spline basis matrix.
% knots1 : array    - Knot sequence of dimension 1.
% knots2 : array    - Knot sequence of dimension 2.
arguments
    X (:,2) double
    nr_splines (1,2) double = [10,10];
    ll (1,2) double = [3,3];
    knot_type (1,2) string = ["e", "e"];
end

    [B1, knots1] = Bspline.basismatrix(X(:,1), nr_splines(1), ll(1), knot_type(1));
    [B2, knots2] = Bspline.basismatrix(X(:,2), nr_splines(2), ll(2), knot_type(2));
    
    T = zeros(size(X,1), prod(nr_splines));
    for i=1:size(X,1)
        T(i,:) = kron(B2(i,:), B1(i,:));
    end
    T = sparse(T);
end

