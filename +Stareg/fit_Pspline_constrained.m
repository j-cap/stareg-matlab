function [coef, B, knots] = fit_Pspline_constrained(X, y, constraint, lam_c, nr_splines, ll, knot_type)
%%
% Iterative 1-D constrained P-spline fit.
%
% Parameters:
% ----------
% X : array           - Input data of shape (n_samples, 1) to calculate the 
%                       coefficents for.
% y : array           - Output data of shape (n_samples, 1)to calculate the 
%                       coefficients for.
% constraint : string - Type of constraint.
% lam_c : double      - Value of the constraint parameter lam_c.
% nr_splines : int    - Number of parameters (== number of B-spline basis 
%                       functions).
% ll : int            - Specifies the order of the B-spline basis functions.
% knot_type : str     - Decide between equidistant "e" and quantile based 
%                       "q" knot placement.
%
% Returns:
% --------
% coef : array   - Least squares coefficients.
% B : matrix     - B-spline basis matrix.
% knots : array  - Knot sequence

arguments
    X (:,1) double
    y (:,1) double
    constraint (1,1) string
    lam_c (1,1) double = 6000
    nr_splines (1,1) double = 10
    ll (1,1) double = 3
    knot_type (1,1) string = "e"
end

    % calc optimal smoothing parameter
    lam_opt = Bspline.calc_GCV(X,y,100, 0, nr_splines, ll, knot_type);
    [coef, B, knots] = Bspline.fit_Pspline(X, y, lam_opt, nr_splines, ll, knot_type);

    % generate mapping matrices
    Ds = Utils.mapping_matrix("smooth", nr_splines);
    Dc = Utils.mapping_matrix(constraint, nr_splines);
    v = Utils.check_constraint(coef, constraint, B, y);
    v_old = zeros(size(v));
    i = 0;
    % iterate till no change in v
    BtB = B'*B;
    DstDs = Ds'*Ds;
    By = B' * y;

    figure()
    scatter(X, y, 'ro', 'DisplayName', 'Data'); hold on;
    plot(X, B*coef, 'DisplayName', 'PLS-Fit', 'LineWidth',2);
    disp(['PLS-Fit MSE = ', num2str(Utils.mse(y, B*coef))]);


    while ~all(v_old == v) 
        v_old = v;
        DVD = Dc'*diag(v)*Dc;
        coef = (BtB + lam_opt * DstDs + lam_c * DVD) \ By;
        v = Utils.check_constraint(coef, constraint, B, y);
        i = i + 1;
        disp(['Iteration: ', num2str(i)]);
        disp(['MSE = ', num2str(Utils.mse(y, B*coef))]);
        plot(X, B * coef, 'DisplayName', ['Iteration: ' num2str(i)], 'LineWidth', 2); hold on
    end
    legend(); grid()
end


