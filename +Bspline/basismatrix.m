function [B, knots] = basismatrix(X, nr_splines, ll, knot_type)
%% 
% Generate the B-spline basis matrix for nr_splines given the data X. 
% 
% Note: (nr_splines + order + 1) knots are needed for a B-spline basis of 
%       order l with nr_splines, 
%       e.g. for l=3, nr_splines=10 -> length(knots) == 14.    
%
% Parameters:
% ----------
% X : array        -  Input data of shape (n_samples, ) to compute the 
%                     B-spline basis matrix for.
% nr_splines : int -  Number of parameters (== number of B-spline basis 
%                     functions).
% ll : int         -  Specifies the order of the B-spline basis functions.
% knot_type : str  -  Decide between equidistant "e" and quantile based "q"
%                     knot placement.
%
% Returns:
% --------
% B  :  matrix   - B-spline basis of dim (length(X) times nsplines).
% knots : array  - Knot sequence.

arguments
   X (:,1) double
   nr_splines (1,1) double
   ll (1,1) double = 3
   knot_type (1,1) string = "e"
end
    
    B = zeros(length(X), nr_splines);
    xmin = min(X); xmax =  max(X);

    if knot_type == "e"
        knots = linspace(xmin, xmax, nr_splines-ll+1);
    elseif knot_type == "q"
        p = linspace(0, 1, nr_splines-ll+1);
        xs = sort(X); % Ordered elements
        i = (length(xs)-1)*p(:)+1; % Indexes
        knots = xs(floor(i))';
    else
        disp("Knot Type not implemented ");
    end
    
    dknots = mean(diff(knots));
    knots = [linspace(xmin-ll*dknots, xmin-dknots, ll), knots, linspace(xmax+dknots, xmax+ll*dknots, ll)];

    for i=ll+1:length(knots)-1
        B(:,i-ll) = Bspline.basisfunction(X, knots, i, ll);
    end
    B = sparse(B);
end