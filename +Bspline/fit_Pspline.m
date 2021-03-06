function [coef, B, knots] = fit_Pspline(X,y, lam, nr_splines, order, knot_type)
%% 
% Penalized Least Squares for B-Splines or tensor-product B-splines.
%
% Parameters:
% ----------
% X : array        - Input data of shape (n_samples, 1 or 2) to calculate the 
%                    coefficents for.
% y : array        - Target data of shape (n_samples, 1)to calculate the 
%                    coefficients for.
% lam : double     - Value of the smoothing parameter lambda.
% nr_splines : int - Number of parameters (== number of B-spline basis 
%                    functions).
% order : int      - Specifies the order of the B-spline basis functions.
% knot_type : str  - Decide between equidistant "e" and quantile-based "q"
%                    knot placement.
%
% Returns:
% --------
% coef : array   - Penalize Least-Squares coefficients.
% B : matrix     - B-spline basis matrix or tensor-product B-spline basis matrix 
% knots : struct - Knot sequences (for both dimensions)
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
   lam (:,1) double = 1;
   nr_splines (:,1) double = 10;
   order (:,1) double = 3;
   knot_type (:,1) string = "e";
end

    [~,c] = size(X);
    if c == 1 % check if B-spline
        [B, knots] = Bspline.basismatrix(X, nr_splines, order, knot_type); % create B-spline basis matrix
        D = Utils.mapping_matrix("smooth", nr_splines); % create mapping matrix
        coef = ((B'*B) + lam * (D'*D)) \ (B' * y);
    elseif c == 2 % or tensor-product B-spline
        [B, knots] = Bspline.tensorproduct_basismatrix(X, nr_splines, order, knot_type);  % create tensor-product B-spline basis matrix
        D1 = Utils.mapping_matrix_tp("smooth", nr_splines, 1); % create mapping matrix for dimension 1
        D2 = Utils.mapping_matrix_tp("smooth", nr_splines, 2); % create mapping matrix for dimension 2
        coef = ((B'*B) + lam * (D1'*D1) + lam * (D2'*D2)) \ (B' * y);
    else
        disp("Maximal dimension  == 2!");
        return
    end
    
end
