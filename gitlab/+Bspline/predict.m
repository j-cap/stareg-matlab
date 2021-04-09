function s = predict(Xpred,knots,coef,order)
%%
% Evaluate the B-spline given the knot sequence in knots with the
% coefficients in coef for data X.
%
% Parameters:
% X : array      - Data points of shape (n_samples, :) to evaluate the 
%                  B-spline on.
% knots : array  - Knot sequence, length = (length(coef) + order + 1).
% coef : array   - B-spline coefficients.
% order : struct - Order of the B-spline basis functions.
%
% Returns
% -------
% s : array      - B-spline values at Xpred.
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
    Xpred (:,:) double
    knots (:,:) struct
    coef (:,1) double
    order (1,:) struct;
end

    [n_data,dims] = size(Xpred);
    if dims == 1
        B = zeros(n_data, length(coef));
        for spline_idx=order.so1+1:length(knots.k1)-1 % evaluate each B-spline basis functions for Xpred
            B(:,spline_idx-order.so1) = Bspline.basisfunction(Xpred, knots.k1, spline_idx, order.so1);
        end
    elseif dims == 2
        B1 = zeros(n_data, length(knots.k1)-1-order.so1); %basis matrix for direction 1
        B2 = zeros(n_data, length(knots.k2)-1-order.so2); %basis matrix for direction 2
        B = zeros(n_data, length(coef));
        
        for spline_idx=order.so1+1:length(knots.k1)-1 % evaluate the B-spline basis functions for dimension 1
            B1(:,spline_idx-order.so1) = Bspline.basisfunction(Xpred(:,1), knots.k1, spline_idx, order.so1);
        end
        for spline_idx=order.so2+1:length(knots.k2)-1 % evaluate the B-spline basis functions for dimension 2
            B2(:,spline_idx-order.so2) = Bspline.basisfunction(Xpred(:,2), knots.k2, spline_idx, order.so2);
        end
        for data_point=1:n_data % calculate the row-wise Kronecker product between the B-spline basis matrices
            B(data_point,:) = kron(B2(data_point,:), B1(data_point,:));
        end
    end
   
    s = B*coef;
end

