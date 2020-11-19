function [b_PLS_GCV, X, best_lam, k] = GCV_PLS(x, y, nsplines, nlam, plot_, analytic)
%% Calculate the generalized cross validation for the given data (x,y)
%
% o) Currently only for 1D!
% o) Spline degree and knot placement are currently hard coded. 
%
% Parameters: 
% -----------
% x  :  array       - Input data.
% y  :  array       - Output data.
% nsplines : int    - Number of splines to use for the GCV.
% nlam     : int    - Number of smoothing paramters to test.
% plot_    : bool   - Plot fit and smoothing paramter curve.
% analytic :        - Use analytic smoothing matrix or FD smoothing matrix.
%
% Returns:
% --------
% b_PLS_GCV : array  - Optimal coefficients for penalized least squares
%                      fit using generalized cross validation.
% X        : matrix  - Used B-spline basis matrix.
% best_lam : float   - Optimal smoothing parameter.
% k        : array   - Knot sequence.
%%
if isrow(x) x = x'; end
if isrow(y) y = y'; end
lams = logspace(-8,8,nlam);
gcvs = zeros(size(lams));

[ndata, ndim] = size(x);

for i=1:length(lams)
    msg = ['try: lam=', num2str(lams(i))];
    disp(msg);
    if ndim == 1
        [beta_PLS, X, S, k] = PLS(x, y, lams(i), nsplines, analytic);
    elseif ndim == 2
        [beta_PLS, X, S] = PLS_tps(x, y, lams(i), nsplines(1), nsplines(2));
    end
    H = X * (((X'*X) + lams(i) * (S'*S)) \ X');
    ypred = X * beta_PLS;
    gcvs(i) = sum( ((y - ypred) ./ (1 - trace(H) / length(y))).^2);
end
[best_gcv, idx] = min(gcvs);
best_lam = lams(idx);
if plot_
    figure(); title("GCV-Search");
    semilogx(lams, gcvs); hold on; 
    xlabel("Lambdas"); ylabel("GCV-Score"); grid on;
end
disp("Best lambda found = "); best_lam
if ndim == 1
    [b_PLS_GCV, X, S, k] = PLS(x, y, best_lam, nsplines, analytic);
elseif ndim == 2
    [b_PLS_GCV, X, S] = PLS_tps(x, y, best_lam, nsplines(1), nsplines(2));
end

end

