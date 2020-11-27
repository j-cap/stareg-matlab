function [best_lam] = calc_GCV(X, y, nlam, plot_, nr_splines, ll, knot_type)
%% 
% Calculate the generalized cross validation for the given data (x,y) 
% according to Fahrmeir, Regression p. 480.
%
% Note: The function works for B-splines and tensor-product B-splines. When
% using tensor-product B-splines, the arguments nr_splines, ll and
% knot_type need be to of size (1x2). 
%
% Parameters:
% ----------
% X : array           - Input data of shape (n_samples, 1 or 2) to calculate 
%                       the coefficents for.
% y : array           - Output data of shape (n_samples, 1)to calculate the 
%                       coefficients for.
% nr_splines : double - Number of parameters (== number of B-spline basis 
%                       functions).
% ll : int            - Specifies the order of the B-spline basis functions.
% knot_type : str     - Decide between equidistant "e" and quantile based 
%                       "q" knot placement.
% nlam : double       - Specifies the number of lambdas to try for the GCV.
% plot_ : bool        - Specifes whether to plot the GCV curve.
%
% Returns:
% --------
% best_lam : double   - Optimal smoothing parameter.

arguments
   X (:,:) double
   y (:,1) double
   nlam (1,1) double = 100
   plot_ (1,1) double = 0
   nr_splines (1,:) double = 10
   ll (1,:) double = 3
   knot_type (1,:) string = "e"
end

    lams = logspace(-8,8,nlam);
    gcvs = zeros(size(lams));
    n_dim = size(X,2);
    
    for i=1:length(lams)
        msg = ['test: lam=', num2str(lams(i))];
        disp(msg);
        if n_dim == 1
            [coef_PLS, B, ~] = Bspline.fit_Pspline(X, y, lams(i), nr_splines, ll, knot_type);
            D = Utils.mapping_matrix("smooth", nr_splines);
            BtB = B'*B;
            traceH = trace(BtB * inv(BtB + lams(i) * (D'*D)));
        elseif n_dim == 2
            [coef_PLS, B, ~] = Bspline.fit_tp_Pspline(X, y, lams(i), nr_splines, ll, knot_type);
            [Id2D1, D2Id1] = Utils.mapping_matrix_tp("smooth", nr_splines);
            BtB = B'*B;
            traceH = trace(BtB * inv(BtB + lams(i)*(Id2D1'*Id2D1) + lams(i)*(D2Id1'*D2Id1)));
        else
            disp("Only 1-D and 2-D supported!");
            return
        end
        ypred = B * coef_PLS;
        gcvs(i) = sum(( (y - ypred) ./ (1 - (traceH/length(y))) ).^2);
    end
    [~, idx] = min(gcvs);
    best_lam = lams(idx);
    if plot_
        figure(); title("GCV-Search");
        set(gcf,'Position',[10 10 1600 1000])
        semilogx(lams, gcvs, 'b'); hold on;
        scatter(lams(idx), gcvs(idx), 'ro', 'DisplayName', 'Optimal Lambda')
        xlabel("Lambdas"); ylabel("GCV-Score"); grid on; legend();
    end
    disp("Best lambda found = "); best_lam

end

