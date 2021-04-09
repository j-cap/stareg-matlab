function [coef_, B, knots] = fit(X, y, nr_splines, order, knot_type)
%%
% Calculate the least squares coefficients for the B-spline or
% tensor-product B-spline fit. 
%
% Inputs:
% ----------
% X : array        - Input data of shape (n_samples, 1 or 2) to calculate  
%                    the coefficents for.
% y : array        - Target data of shape (n_samples, 1) to calculate the 
%                    coefficients for.
% nr_splines : int - Number of parameters (== number of B-spline basis 
%                    functions).
% order : int      - Specifies the order of the B-spline basis functions.
% knot_type : str  - Decide between equidistant "e" and quantile-based "q"
%                    knot placement.
%
% Outputs:
% --------
% coef : array   - Least squares coefficients.
% B : matrix     - B-spline or tensor-product B-spline basis matrix.
% knots : array  - Knot sequence.
%
% Dependencies:
%    Matlab release: 2020b
%
% This function is part of: stareg-matlab
%
% Author:  Jakob Weber
% email:   jakob.weber@ait.ac.at
% Company: Austrian Institute of Technology GmbH
%          Complex Dynamical Systems
%          Center for Vision, Automation & Control
%          http://www.ait.ac.at
%
% Version: 1.0.1 - 2021-02-02

% Change log:
% x.y.z - 2021-02-02 - author:
% - added important feature, s. issue #34
% - fixed bug #2
%%
arguments % default values
    X (:,:) double
    y (:,1) double
    nr_splines (:,:) double = 10;
    order (:,:) double = 3;
    knot_type (:,:) string = "e";
end
    
    [~, cols] = size(X);
    if cols == 1 % check if B-spline
        [B, knots] = Bspline.basismatrix(X, nr_splines, order, knot_type);
    elseif cols == 2 % or tensor-produc B-spline
        [B, knots] = Bspline.tensorproduct_basismatrix(X, nr_splines, order, knot_type);
    else
        disp('Maximal dimension == 2!');
        return;
    end
    coef_ = B \ y; % Least-Squares 

end

