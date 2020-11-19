function [beta_PLS, X, S] = PLS_tps(x,y,lam, nsplines1, nsplines2)
%%
% Penalized Least Squares for B-Splines -> P-Splines
%
% Parameters:
% -----------
% x   : array  - predictors
% y   : array  - target values
% lam : float  - smoothing parameter
% nsplines1 : int  - Nr. of splines to use for dimension 1
% nsplines2 : int  - Nr. of splines to use for dimension 2
%
%
% Returns:
% --------
% beta_PLS : array   - computed coefficients
% X        : matrix  - tp-spline basis
% S        : matrix  - smoothing matrix
%%

X = tpspline_basis(x, nsplines1, nsplines2, 2);
S = smoothness_matrix_TPS(nsplines1, nsplines2, 2);
beta_PLS = ((X'*X) + lam * (S'*S)) \ X' * y;
end
