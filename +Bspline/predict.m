function s = predict(X,knots,coef,ll)
%%
% Evaluate the B-spline given the knot sequence in knots with the
% coefficients in coef.
%
% Parameters:
% X : array      - Data points of shape (n_samples, 1) to evaluate the 
%                  B-spline on.
% knots : array  - Knot sequence.
% coef : array   - B-spline coefficients.
% ll : double    - Order of the B-spline basis functions.
%
% Returns
% -------
% s : array      - B-spline values at X.

arguments
    X (:,1) double
    knots (:,1) double
    coef (:,1) double
    ll (1,1) double = 3
end

% generate B-spline basis given the knots for the values in X
B = zeros(length(X), length(coef));
for i=ll+1:length(knots)-1 
    B(:,i-ll) = Bspline.basisfunction(X, knots, i, ll);
end
s = B*coef

end

