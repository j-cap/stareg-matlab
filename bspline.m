function [y] = bspline(x, knots, j, l)    
%% B-Splines according to Fahrmeir, Regression p. 429
%
% Parameters:
% -----------
% x : array
%     Data to calculate the spline values for.
% k : array
%     Array of knot locations.
% j : int
%     Index of the b-spline basis function to compute.
% l : int
%     Order of the spline.
%
% Returns:
% --------
% y : B-spline evaluated at x.

    if l==0
        y = knots(j) <= x & x < knots(j+1);
    else
        z0 = (x - knots(j-l)) / (knots(j) - knots(j-l));
        z1 = (knots(j+1) - x) / (knots(j+1) - knots(j+1-l));
        y = z0.*bspline(x, knots, j-1, l-1) + z1.*bspline(x, knots, j, l-1);
    end

end



