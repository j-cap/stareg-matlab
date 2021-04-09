function [B, knots] = basismatrix(X, nr_splines, order, knot_type)
%% 
% Generate the B-spline basis matrix for nr_splines given the data X. 
% 
% Note: (nr_splines + order + 1) knots are needed for a B-spline basis of 
%       'order' with nr_splines, 
%       e.g. for order=3, nr_splines=10 -> length(knots) == 14.    
%
% Inputs:
% ----------
% X : array        -  Input data of shape (n_samples, ) to compute the 
%                     B-spline basis matrix for.
% nr_splines : int -  Number of parameters (== number of B-spline basis 
%                     functions).
% order : int      -  Specifies the order of the B-spline basis functions.
% knot_type : str  -  Decide between equidistant "e" and quantile-based "q"
%                     knot placement.
%
% Outputs:
% --------
% B  :  matrix    - B-spline basis of dim (length(X) x nr_splines).
% knots : struct  - Knot sequence for dimension 1 (.k1).
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
   X (:,1) double;
   nr_splines (1,1) double = 10;
   order (1,1) double = 3;
   knot_type (1,1) string = "e";
end
    
    B = zeros(length(X), nr_splines);
    xmin = min(X);
    xmax = max(X);

    if knot_type == "e"
        inner_knots = linspace(xmin, xmax, nr_splines-order+1);
    elseif knot_type == "q"
        p = linspace(0, 1, nr_splines-order+1);
        xs = sort(X); 
        i = (length(xs)-1)*p(:)+1; 
        inner_knots = xs(floor(i))';
    else
        disp("Knot Type not implemented ");
    end
    
    dknots = mean(diff(inner_knots));
    knots_left = linspace(xmin-order*dknots, xmin-dknots, order);
    knots_right = linspace(xmax+dknots, xmax+order*dknots, order);
    knots = [knots_left, inner_knots, knots_right]; 

    for spline_idx=order+1:length(knots)-1
        B(:,spline_idx-order) = Bspline.basisfunction(X, knots, spline_idx, order);
    end
    
    knots = struct("k1", knots);
end