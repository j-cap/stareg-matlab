function [b_PLS, X, S, k] = PLS(x,y,lam, nsplines, analytic)
%% Penalized Least Squares for B-Splines 
%
% o) Spline order and knot_type are hard coded in line 22.
%
% Parameters:
% -----------
% x   : array  - Input variables.
% y   : array  - Output variables.
% lam : float  - Smoothing parameter.
% nsplines : int  - Number of splines.
% analytic : bool - Specify the use of analytic or FD smoothness matrix.
%
% Returns:
% --------
% b_PLS : array   - computed coefficients
% X        : matrix  - b-spline basis
% S        : matrix  - smoothing matrix
% k        : array   - knot sequence
%%
    if isrow(x) x = x'; end 
    if isrow(y) y = y'; end
    [X, k] = bspline_basis(x, nsplines, 3, "q");
    if analytic == 0
        S = smoothness_matrix(nsplines, 2, 0);
    else
        S = smoothness_matrix(nsplines, 2, k);
    end
    b_PLS = ((X'*X) + lam * (S'*S)) \ (X' * y);
end
