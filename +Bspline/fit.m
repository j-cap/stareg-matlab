function [coef_, B, knots] = fit(X, y, nr_splines, order, knot_type)
%%
% Calculate the least squares coefficients for the B-spline or
% tensor-product B-spline fit. 
%
%
% Parameters:
% ----------
% X : array        - Input data of shape (n_samples, 1 or 2) to calculate  
%                    the coefficents for.
% y : array        - Target data of shape (n_samples, 1)to calculate the 
%                    coefficients for.
% nr_splines : int - Number of parameters (== number of B-spline basis 
%                    functions).
% order : int      - Specifies the order of the B-spline basis functions.
% knot_type : str  - Decide between equidistant "e" and quantile based "q"
%                    knot placement.
%
% Returns:
% --------
% coef : array   - Least squares coefficients.
% B : matrix     - B-spline or tensor-product B-spline basis matrix.
% knots : array  - Knot sequence.

arguments
    % default values are for B-spline fits
    X (:,:) double
    y (:,1) double
    nr_splines (:,:) double = 10;
    order (:,:) double = 3;
    knot_type (:,:) string = "e";
end
    
    [~, cols] = size(X);
    if cols == 1
        [B, knots] = Bspline.basismatrix(X, nr_splines, order, knot_type);
    elseif cols == 2
        [B, knots] = Bspline.tensorproduct_basismatrix(X, nr_splines, order, knot_type);
    else
        disp('Maximal dimension == 2!');
        return;
    end
    coef_ = B \ y;

end

