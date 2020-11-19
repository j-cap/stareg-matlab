function [b, X_pls] = PIRLS(x, y, nsplines, constraint, plot_, analytic)
%% 1D Penalized Iteratively Reweighted Least Squares according to Hofner 2011.
%
%
% Parameters:
% -----------
% x  : array - input data
% y  : array - output data
% nsplines : int - Nr. of splines to use
% constraint : string - Type of the constraint
% plot_      : bool  - Show plots
% analytic   : bool  - Use analytic derivative or FD smoothness matrix
%%
if isrow(x) x = x'; end
if isrow(y) y = y'; end 

if plot_
    figure(); hold on; scatter(x, y); legend();
end

% calc optimal smoothing parameter
[b, X_pls, lam_pls, k] = GCV_PLS(x,y,nsplines,100,0, analytic);

if plot_
    plot(x, X_pls * b, 'DisplayName', 'PLS Estimate', 'LineWidth', 2);
end

% generate mapping matrices
if analytic == 0
    use_knots = 0;
else
    use_knots = k;
end
S = smoothness_matrix(nsplines, 2, use_knots);
D = mapping_matrix(constraint, nsplines, use_knots);
v = check_constraint(b, constraint, y, X_pls);
v_old = zeros(size(v));
i = 0;
% iterate till no change in v
while ~all(v_old == v)
    i = i + 1;
    disp(['Iteration: ', num2str(i)]);
    v_old = v;
    XTX = X_pls'*X_pls;
    STS = S'*S;
    DVD = D'*diag(v)*D;
    Xy = X_pls' * y;
    b = (XTX + lam_pls * STS + 6000 * DVD) \ Xy;
    v = check_constraint(b, constraint, y, X_pls);
    plot(x, X_pls * b, 'DisplayName', ['Iteration: ' num2str(i)], 'LineWidth', 2); hold on
end


