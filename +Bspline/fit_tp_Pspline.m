function [coef, T, knots1, knots2] = fit_tp_Pspline(X,y, lam, nr_splines, ll, knot_type)
%% 
% Penalized Least Squares for B-Splines 
%
% Parameters:
% ----------
% X : array          - Input data of shape (n_samples, 2) to calculate the 
%                      coefficents for.
% y : array          - Output data of shape (n_samples, 1)to calculate the 
%                      coefficients for.
% lam : double       - Value of the smoothing parameters lambda.
% nr_splines : array - Number of parameters (== number of B-spline basis 
%                      functions).
% ll : array         - Specifies the orders of the B-spline basis functions.
% knot_type : array  - Decide between equidistant "e" and quantile based "q"
%                      knot placement.
%
% Returns:
% --------
% coef : array    - Least squares coefficients.
% T : matrix      - B-spline basis matrix.
% knots1 : array  - Knot sequence for dimension 1.
% knots2 : array  - Knot sequence for dimension 2.

arguments
   X (:,2) double
   y (:,1) double
   lam (1,1) double = 1
   nr_splines (1,2) double = [10,10]
   ll (1,2) double = [3,3]
   knot_type (1,2) string = ["e","e"]
end
    [T, knots1, knots2] = Bspline.tp_basismatrix(X, nr_splines, ll, knot_type);
    [Id2D1, D2Id1] = Utils.mapping_matrix_tp("smooth", nr_splines);
    coef = ((T'*T) + lam * (Id2D1'*Id2D1) + lam * (D2Id1'*D2Id1)) \ (T' * y);
end
