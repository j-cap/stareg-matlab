function [T,knots] = tensorproduct_basismatrix(X, nr_splines, orders, knot_types)
%% 
% Generate the tensor-product B-spline basis given the data X!
% 
% Parameters
% ----------
% X : matrix            - Input data of shape (n_samples, 2).
% nr_splines : array    - Number of parameters for dimension 1 and 2.
% orders : array        - Order of the B-spline basis functions for 
%                         each dimension.
% knot_types : array    - Knot type of the B-spline basis functions 
%                         for each dimension 
%
% Returns:
% --------
% T  : matrix       - Tensor-product B-spline basis matrix.
% knots : struct    - Knot sequence for dimension 1 (.k1) and 2 (.k2).

arguments
    X (:,2) double
    nr_splines (1,2) double = [10,10];
    orders (1,2) double = [3,3];
    knot_types (1,2) string = ["e", "e"];
end

    [B1, knots1] = Bspline.basismatrix(X(:,1), nr_splines(1), orders(1), knot_types(1));
    [B2, knots2] = Bspline.basismatrix(X(:,2), nr_splines(2), orders(2), knot_types(2));
    
    T = zeros(size(X,1), prod(nr_splines));
    for i=1:size(X,1)
        T(i,:) = kron(B2(i,:), B1(i,:));
    end
    
    knots.k1 = knots1.k1;
    knots.k2 = knots2.k1;
end

