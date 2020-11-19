function [X, knots] = bspline_basis(x, nsplines, sorder, knot_type)
%% Generate the B-spline basis for nsplines given the data x!
% 
% o) One needs nsplines + order + 1 knots for a spline basis of order m with k parameters. 
% o) Equidistant and quantile based knot placemnet are available.
%
%
% Parameters:
% ----------
% x : np.ndarray   -  Data of shape (n_samples, ) to compute the B-spline with.
% nsplines : int   -  Number of parameters (== number of B-splines).
% sorder : int     -  Specifies the order of the spline.
% knot_type : bool -  Decide betweend equidistant "e" and quantile based "q"
%                     knot placement.
%
% Returns:
% --------
% X  :  matrix   - B-spline basis of dim length(x) times nsplines.
% knots : array  - Knot sequence.
%%
    ndata = length(x);
    X = zeros(ndata, nsplines);

    xmin = min(x);
    xmax =  max(x);

    if knot_type == "e"
        knots = linspace(xmin, xmax, nsplines-sorder+1);
    elseif knot_type == "q"
        p = linspace(0, 1, nsplines-sorder+1);
        xs = sort(x); % Ordered elements
        [n,m] = size(xs);
        i = (n-1)*p(:)+1; % Indexes
        knots = xs(floor(i))';
    else
        disp("Knot Type not implemented ");
    end
    
    dknots = mean(diff(knots));
    knots = [linspace(xmin-sorder*dknots, xmin-dknots, sorder), knots, linspace(xmax+dknots, xmax+sorder*dknots, sorder)];

    for i=sorder+1:length(knots)-1
        X(:,i-sorder) = bspline(x, knots, i, sorder);
    end
    X = sparse(X);
end