function [T,knots] = tensorproduct_basismatrix(X, nr_splines, orders, knot_types)
%% 
% Generate the tensor-product B-spline basis given the data X!
% 
% Inputs:
% ----------
% X : matrix            - Input data of shape (n_samples, 2).
% nr_splines : array    - Number of parameters for dimension 1 and 2.
% orders : array        - Order of the B-spline basis functions for 
%                         each dimension.
% knot_types : array    - Knot type of the B-spline basis functions 
%                         for each dimension 
%
% Outputs:
% --------
% T  : matrix       - Tensor-product B-spline basis matrix.
% knots : struct    - Knot sequence for dimension 1 (.k1) and 2 (.k2).
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
    X (:,2) double
    nr_splines (1,2) double = [10,10];
    orders (1,2) double = [3,3];
    knot_types (1,2) string = ["e", "e"];
end

    % create the B-spline basis matrices in both dimensions
    [B1, knots1] = Bspline.basismatrix(X(:,1), nr_splines(1), orders(1), knot_types(1));
    [B2, knots2] = Bspline.basismatrix(X(:,2), nr_splines(2), orders(2), knot_types(2));
    
    [n_data,~] = size(X);
    T = zeros(n_data, prod(nr_splines));
    for data_point=1:n_data % calculate the row-wise Kronecker product between the B-spline basis matrices
        T(data_point,:) = kron(B2(data_point,:), B1(data_point,:));
    end
    
    knots.k1 = knots1.k1;
    knots.k2 = knots2.k1;
end

