function b = basisfunction(X, knots, j, ll)    
%% 
% B-Spline basis function definition according to 
% Fahrmeir, Regression p. 429
%
% Parameters:
% -----------
% X :  array      - Input data of shape (n_samples, 1) to evaluate the 
%                   B-spline basis function of order l on.
% knots :  array  - Knot sequence defining the B-spline basis function.
% j :  int        - Index of the B-spline basis function to evaluate.
% ll : int        - Order of the B-spline basis function, e.g. l=3 -> cubic
%
% Returns:
% --------
% b : array    - B-spline evaluated at x.

arguments
    X (:,1) double 
    knots (:,1) double
    j (1,1) double
    ll (1,1) double = 3
end
    if ll==0
        b = knots(j) <= X & X < knots(j+1);
    else
        z0 = (X - knots(j-ll)) / (knots(j) - knots(j-ll));
        z1 = (knots(j+1) - X) / (knots(j+1) - knots(j+1-ll));
        b = z0.*Bspline.basisfunction(X, knots, j-1, ll-1) + z1.*Bspline.basisfunction(X, knots, j, ll-1);
    end

end
