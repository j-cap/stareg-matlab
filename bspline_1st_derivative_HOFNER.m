function [dSpline, dX] = bspline_1st_derivative_HOFNER(x, knots, coef, l)
%% Calculate the first derivative according to Fahrmeir, Regression p. 430
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
 
    dX = zeros(length(x), length(coef)-1);
    for i=l+1:length(knots)-2
        dX(:,i-l) = bspline(x, knots, i, l-1);
    end
    h = mean(diff(knots));
    dcoef = diff(coef);
    dSpline = (1/h) * (dX * dcoef);
    
end

