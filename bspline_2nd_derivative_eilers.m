function [dSpline, dX] = bspline_2nd_derivative_eilers(x, knots, coef, l)
%% Calculate the first derivative according to Eilers and Marx, 1995
%       !!! for equidistant grids !!!
%
% Parameters:
% ----------
% x : array     - input data vector
% knots : array     - knot sequence for the B-splines
% coef  : array     - estimated coefficients of the B-splines
% l     : int       - order of the B-splines
%
% Returns:
% --------
% dSpline  : array       - Values of the derivative at x
% dX       : matrix      - derivative of the splines
%%
 
    dX = zeros(length(x), length(coef)-2);
    for i=1:length(coef)-2
        dX(:,i) = bspline(x, knots, i+l, l-2);
    end
    dcoef = diff(coef, 2);
    h = mean(diff(knots));
    dSpline = (dX * dcoef) / h^2;
    
    dX = h^2;
end

