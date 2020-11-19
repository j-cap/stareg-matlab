function [X,knots1,knots2] = tpspline_basis(xdata,nsplines1,nsplines2,sorder)
%% Generate the Tensorproduct-spline basis given the data xdata!
% 
% o) The knot placement is hard coded in lines 16 and 17.
%
% Parameters
% ----------
% nsplines1 : int
%     Number of parameters for dimension 1(== number of B-splines).
% nsplines2 : int
%     Number of parameters for dimension 2(== number of B-splines).
% sorder : int, "spline order"
%     Specifies the order of the spline.
% xdata : np.ndarray
%     Data of shape (n_samples, 2) to compute the B-spline with.
%
% Returns:
% --------
% X  : matrix       - Tensorproduct-Spline basis.
% knots1 : array    - Knot sequence of dimension 1.
% knots2 : array    - Knot sequence of dimension 2.
%%
    [B1, knots1] = bspline_basis(xdata(:,1), nsplines1, sorder, "e");
    [B2, knots2] = bspline_basis(xdata(:,2), nsplines2, sorder, "e");
    
    n = size(xdata);
    X = zeros(n(1), nsplines1*nsplines2);
    for i=1:n(1)
        X(i,:) = kron(B2(i,:), B1(i,:));
    end
    X = sparse(X);
end

