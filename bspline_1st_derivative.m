function [dSpline, dX] = bspline_1st_derivative(x, knots, coef, l)
%% Calculate the first derivative according to Fahrmeir, Regression p. 430
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
    h = zeros(length(coef)-1,1);
    %for i=l+1:length(knots)-2
    for i=2:length(coef)
        dX(:,i-1) = Bspline.basisfunction(x, knots, l+i-1, l-1);
        h(i-1) = knots(l+i) - knots(i); 
    end
    dcoef = diff(coef) ./ h;
    dSpline = l * (dX * dcoef);
    
end

