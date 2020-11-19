function S = smoothness_matrix(nparam, degree, knots)
%% Calculate the smoothness matrix for PIRLS and PLS. 
%
% We have currently 2 options determined by the parameter knots: 
%  1. if knots == 0, we use the FD-approximation given in Fahrmeir,
%     Regression p.438ff
%  2. if knots == knot sequence, we use the real derivative values in the
%     smoothness matrix, 
%
% Parameters:
% ------------
%   nparam : int  - Number of splines. 
%   degree : int  - Degree of the smoothness matrix, either 1 or 2.
%   knots  : array or 0   - Knot vector or 0.
%
% Returns:
% ------------
%   S     : matrix - Smoothness penalty matrix.
%
% o) For FD-approx, S is of shape nparam-degree times nparam
% o) For real derivative, S is of shape nparam times length(knots) - In
% each column, the values for the 2nd derivative are given for this knot.
%%

    if knots == 0
        e = ones(nparam, 1);
        if degree == 1
            diagonals = [-e, e];
        elseif degree == 2
            diagonals = [e, -2*e, e];
        end
        S = spdiags(diagonals, 0:degree, nparam-degree, nparam);
    else
        l = length(knots) - nparam - 1;
        dS = zeros(nparam, length(knots));
        % calculate first derivative according to Fahrmeir, Regression
        % p.429
        if degree == 1
            for i=l+1:length(knots)-1
                left = bspline(knots, knots, i-1, l-1) / (knots(i) - knots(i-l)); 
                right = bspline(knots, knots, i, l-1) / (knots(i+1) - knots(i+1-l));
                dS(i-l,:) = l * (left - right);
            end 
            S = sparse(dS');
        elseif degree == 2
            for i=l+1:length(knots)-1
                left = bspline(knots, knots, i-1, l-1) / (knots(i) - knots(i-l));
                right = bspline(knots, knots, i, l-1) / (knots(i+1) - knots(i+1-l));
                dS(i-l,:) = l * (left - right);
            end 
            % calculate second derivative using FD
            ddS = diff(dS')./diff(knots');
            S = sparse(ddS);
        end
    end
end
