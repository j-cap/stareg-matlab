function [dSpline, dX, dcoef, h2] = bspline_2nd_derivative_fahrmeir(x, knots, coef, l)
%% Calculate the first derivative according to Fahrmeir, Regression
%       !!! for both grid types grids !!!
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
    h2 = zeros(length(coef)-2, 1);
    dcoef = zeros(length(coef)-2,1);
    for i=3:length(coef)
        dX(:,i-2) = Bspline.basisfunction(x, knots, i+l-2, l-2);
        h2(i-2) = (knots(i+l) - knots(i)) * (knots(i+l+1) - knots(i+1));
        dcoef(i-2) = (coef(i) - coef(i-1)) * (coef(i-1) - coef(i-2));
    end
    dcoef = diff(coef, 2);
    dSpline =  l^2 * dX * (dcoef ./ h2);
    
end

