function [coef, B, knots] = fit_Pspline(X,y, lam, nr_splines, ll, knot_type)
%% 
% Penalized Least Squares for B-Splines 
%
% Parameters:
% ----------
% X : array        - Input data of shape (n_samples, 1) to calculate the 
%                    coefficents for.
% y : array        - Output data of shape (n_samples, 1)to calculate the 
%                    coefficients for.
% lam : double     - Value of the smoothing parameter lambda.
% nr_splines : int - Number of parameters (== number of B-spline basis 
%                    functions).
% ll : int         - Specifies the order of the B-spline basis functions.
% knot_type : str  - Decide between equidistant "e" and quantile based "q"
%                    knot placement.
%
% Returns:
% --------
% coef : array   - Least squares coefficients.
% B : matrix     - B-spline basis matrix.
% knots : array  - Knot sequence.
arguments
   X (:,:) double
   y (:,1) double
   lam (:,1) double = 1;
   nr_splines (:,1) double = 10;
   ll (:,1) double = 3;
   knot_type (:,1) string = "e";
end

    [~,c] = size(X);
    if c == 1
        [B, knots] = Bspline.basismatrix(X, nr_splines, ll, knot_type);
        D = Utils.mapping_matrix("smooth", nr_splines);
        coef = ((B'*B) + lam * (D'*D)) \ (B' * y);
    elseif c == 2
        [B, knots1, knots2] = Bspline.tensorproduct_basismatrix(X, nr_splines, ll, knot_type);
        knots = {knots1, knots2};
        D1 = Utils.mapping_matrix_tp("smooth", nr_splines, 1);
        D2 = Utils.mapping_matrix_tp("smooth", nr_splines, 2);
        coef = ((B'*B) + lam * (D1'*D1) + lam * (D2'*D2)) \ (B' * y);
    else
        disp("Maximal dimension  == 2!");
        return
    end
    
end
