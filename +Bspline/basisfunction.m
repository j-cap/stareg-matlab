function b = basisfunction(X, knots, knot_idx, order)    
%% 
% B-Spline basis function definition according to 
% Fahrmeir, Regression p. 429
%
% Parameters:
% -----------
% X :  array        - Input data of shape (n_samples, 1) to evaluate the 
%                     B-spline basis function of 'order' on.
% knots :  array    - Knot sequence defining the B-spline basis function.
% knot_idx :  int   - Index of the B-spline basis function to evaluate.
% order : int       - Order of the B-spline basis function, e.g. 3 -> cubic
%
% Returns:
% --------
% b : array    - B-spline evaluated at x.

arguments
    X (:,1) double 
    knots (:,1) double
    knot_idx (1,1) double
    order (1,1) double = 3
end
    if order==0
        b = knots(knot_idx) <= X & X < knots(knot_idx+1);
    else
        z0 = (X - knots(knot_idx-order)) / (knots(knot_idx) - knots(knot_idx-order));
        z1 = (knots(knot_idx+1) - X) / (knots(knot_idx+1) - knots(knot_idx+1-order));
        b1 = z0.*Bspline.basisfunction(X, knots, knot_idx-1, order-1);
        b2 = z1.*Bspline.basisfunction(X, knots, knot_idx, order-1);
        b = b1 + b2;
    end

end
