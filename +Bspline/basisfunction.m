function b = basisfunction(X, knots, spline_idx, order)    
%% 
% B-Spline basis function definition according to Fahrmeir, Regression p. 429
%
% Inputs:
% -----------
% X :  array        - Input data of shape (n_samples, 1) to evaluate the 
%                     B-spline basis function of 'order' on.
% knots :  array    - Knot sequence defining the B-spline basis function.
% spline_idx :  int   - Index of the B-spline basis function to evaluate.
% order : int       - Order of the B-spline basis function, e.g. 3 -> cubic.
%
% Outputs:
% --------
% b : array    - B-spline basis function evaluated at x.
%
% % Dependencies:
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

arguments
    X (:,1) double 
    knots (:,1) double
    spline_idx (1,1) double
    order (1,1) double = 3
end
    if order==0
        b = knots(spline_idx) <= X & X < knots(spline_idx+1);
    else
        z0 = (X - knots(spline_idx-order)) / (knots(spline_idx) - knots(spline_idx-order));
        z1 = (knots(spline_idx+1) - X) / (knots(spline_idx+1) - knots(spline_idx+1-order));
        b1 = z0.*Bspline.basisfunction(X, knots, spline_idx-1, order-1);
        b2 = z1.*Bspline.basisfunction(X, knots, spline_idx, order-1);
        b = b1 + b2;
    end

end
