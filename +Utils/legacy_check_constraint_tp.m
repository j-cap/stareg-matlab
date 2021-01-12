function v = check_constraint(coef, constraint, nr_splines, dim)
%% 
% Check whether the coefficients in coef hold the constraint
%
% Parameters:
% -----------
% coef  : array       - Array of coefficients.
% constraint : str    - Constraint type.
% nr_splines : array  - Nr of basis functions for each dimension.
% dim : double        - Indicator variable for the dimension of the
%                       constraint to check the coef against.
%
% Returns:
% --------
% v  : array      - Diagonal elements of the weighting matrix V.

arguments
    coef (:,1) double
    constraint (1,1) string = "inc"
    nr_splines (1,2) double = [10,10]
    dim (1,1) double = 1
end

    threshold = 0;
    coef = reshape(coef, nr_splines(1), nr_splines(2));
    if constraint == "inc"
        v = diff(coef, 1, dim) < threshold;
    elseif constraint == "dec"
        v = diff(coef, 1, dim) > -threshold;
    elseif constraint == "conv"
        v = diff(coef, 2, dim) < threshold;
    elseif constraint == "conc"
        v = diff(coef, 2, dim) > -threshold;
    elseif constraint == "none"
        v = zeros(size(coef));
    else
        disp("Constraint of type -"+constraint+"- not implemented");
    end
    
    v = v(:);
end
