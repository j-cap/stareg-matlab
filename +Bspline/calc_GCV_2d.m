function [best_lam] = calc_GCV_2d(X, y, nlam, plot_, nr_splines, order, knot_type)
%% 
% Calculate the generalized cross validation for the given data (x,y) 
% according to Fahrmeir, Regression p. 480.
%
% Note: The function works for B-splines and tensor-product B-splines. When
% using tensor-product B-splines, the arguments nr_splines, order and
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
   X (:,2) double
   y (:,1) double
   nlam (1,1) double = 100
   plot_ (1,1) double = 0
   nr_splines (1,2) double = 10
   order (1,2) double = [3,3];
   knot_type (1,2) string = ["e", "e"];
end

    lams = logspace(-8,8,nlam);
    gcvs = zeros(size(lams));
    
    [basismatrix, knots] = Bspline.tensorproduct_basismatrix(X, nr_splines, order, knot_type);
    D1 = Utils.mapping_matrix_tp("smooth", nr_splines, 1);
    D2 = Utils.mapping_matrix_tp("smooth", nr_splines, 2);
    
    BtB = basismatrix'*basismatrix;
    Bty = basismatrix'*y;
    D1tD1 = D1'*D1;
    D2tD2 = D2'*D2;
    
    for i=1:length(lams)
        msg = ['try: lam=', num2str(lams(i))];
        disp(msg);
        coef_pls = (BtB + lams(i) * D1tD1 + lams(i) * D2tD2) \ Bty;
        traceH = trace((BtB + lams(i) * D1tD1 + lams(i) * D2tD2) \ BtB);
        ypred = basismatrix * coef_pls;
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

