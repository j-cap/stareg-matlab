function [best_lam] = calc_GCV(X, y, nlam, plot_, nr_splines, order, knot_type)
%% 
% Calculate the generalized cross validation for the given data (x,y) 
% according to Fahrmeir, Regression p. 480.
%
%
% Parameters:
% ----------
% X : array           - Input data of shape (n_samples, 1) to calculate 
%                       the coefficents for.
% y : array           - Target data of shape (n_samples, 1)to calculate the 
%                       coefficients for.
% nr_splines : double - Number of parameters (== number of B-spline basis 
%                       functions).
% order : int         - Specifies the order of the B-spline basis functions.
% knot_type : str     - Decide between equidistant "e" and quantile based 
%                       "q" knot placement.
% nlam : double       - Specifies the number of lambdas to try for the GCV.
% plot_ : bool        - Specifes whether to plot the GCV curve.
%
% Returns:
% --------
% best_lam : double   - Optimal smoothing parameter.

arguments
   X (:,1) double
   y (:,1) double
   nlam (1,1) double = 100;
   plot_ (1,1) double = 0;
   nr_splines (1,1) double = 10;
   order (1,1) double = 3;
   knot_type (1,1) string = "e";
end

    lams = logspace(-8,8,nlam);
    gcvs = zeros(size(lams));
    
    [B, knots] = Bspline.basismatrix(X, nr_splines, order, knot_type);
    D = Utils.mapping_matrix("smooth", nr_splines);
    
    BtB = B'*B;
    Bty = B'*y;
    DtD = D'*D;
    
    for i=1:length(lams)
        msg = ['try: lam=', num2str(lams(i))];
        disp(msg);
        coef_pls = (BtB + lams(i) * DtD) \ Bty;
        traceH = trace((BtB + lams(i) * DtD) \ BtB);
        ypred = B * coef_pls;
        gcvs(i) = sum(((y - ypred) ./ (1 - traceH/length(y))).^2) / length(y);
    end
    [~, idx] = min(gcvs);
    best_lam = lams(idx);
    if plot_
        figure(); title("GCV-Search");
        set(gcf,'Position',[10 10 1600 1000]);
        semilogx(lams, gcvs, 'b'); hold on;
        scatter(lams(idx), gcvs(idx), 'ro', 'DisplayName', 'Optimal Lambda');
        xlabel("Lambdas"); ylabel("GCV-Score"); grid on; legend();
    end
    disp("Best lambda found = "); best_lam

end

