function [coef, B, knots1, knots2] = fit_tp_Pspline_constrained(X, y, constraint, lam_c, nr_splines, ll, knot_type)
%%
% Iterative 1-D constrained P-spline fit.
%
% Parameters:
% ----------
% X : array           - Input data of shape (n_samples, 2) to calculate the 
%                       coefficents for.
% y : array           - Output data of shape (n_samples, 1)to calculate the 
%                       coefficients for.
% constraint : string - Type of constraint.
% lam_c : double      - Value of the constraint parameter lam_c.
% nr_splines : array  - Number of parameters (== number of B-spline basis 
%                       functions).
% ll : arry           - Specifies the order of the B-spline basis functions.
% knot_type : array   - Decide between equidistant "e" and quantile based 
%                       "q" knot placement.
%
% Returns:
% --------
% coef : array   - Least squares coefficients.
% B : matrix     - B-spline basis matrix.
% knots1 : array  - Knot sequence for dimension 1.
% knots2 : array  - Knot sequence for dimension 2.

arguments
    X (:,2) double
    y (:,1) double
    constraint (1,1) string
    lam_c (1,1) double = 6000
    nr_splines (1,2) double = [10,10]
    ll (1,2) double = [3,3]
    knot_type (1,2) string = ["e","e"]
end

    % calc optimal smoothing parameter
    lam_opt = Bspline.calc_GCV(X,y,100, 0, nr_splines, ll, knot_type);
    [coef, B, knots1, knots2] = Bspline.fit_tp_Pspline(X, y, lam_opt, nr_splines, ll, knot_type);
  
    % generate mapping matrices
    [Ds1, Ds2] = Utils.mapping_matrix_tp("smooth", nr_splines);
    [Dc1, Dc2] = Utils.mapping_matrix_tp(constraint, nr_splines);
    v1 = Utils.check_constraint_tp(coef, constraint, nr_splines, 1);
    v2 = Utils.check_constraint_tp(coef, constraint, nr_splines, 2);
    v = [v1, v2];
    v_old = zeros(size(v));
    i = 0;
    % iterate till no change in v
    BtB = B'*B;
    DstDs_1 = Ds1'*Ds1;
    DstDs_2 = Ds2'*Ds2;    
    By = B' * y;

    disp(['PLS-Fit MSE = ', num2str(Utils.mse(y, B*coef))]);
    while ~all(v_old == v) 
        v_old = v;
        DVD_1 = Dc1'*diag(v1)*Dc1;
        DVD_2 = Dc2'*diag(v2)*Dc2;

        coef = (BtB + lam_opt * DstDs_1 + lam_opt * DstDs_2 + lam_c * DVD_1 + lam_c * DVD_2) \ By;
        v1 = Utils.check_constraint_tp(coef, constraint, nr_splines, 1);
        v2 = Utils.check_constraint_tp(coef, constraint, nr_splines, 2);
        v = [v1, v2];
        i= i + 1;
        disp(['Iteration: ', num2str(i)]);
        disp(['MSE = ', num2str(Utils.mse(y, B*coef))]);
    end
end


