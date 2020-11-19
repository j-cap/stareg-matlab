function v = check_constraint(coef, constraint, y, X)
%% Check whether the coefficients in coef hold the constraint
%
% o) Currently only for 1D!
%
% Parameters:
% -----------
% coef  : array     - Array of coefficients.
% constraint : str  - Constraint type.
% y  : array        - Output data.
% X  : matrix       - B-spline basis.
%
% Returns:
% --------
% v  : array      - Diagonal elements of the weighting matrix V.
%
%%
    if isrow(coef), coef = coef'; end

    threshold = 1e-8;
    if constraint == "inc"
        disp("Check Increasing Constraint");
        v = diff(coef) < threshold;
    elseif constraint == "dec"
        disp("Check Decreasing Constraint");
        v = diff(coef) > -threshold;
    elseif constraint == "conv"
        disp("Check Convexe Constraint");
        v = diff(coef, 2) < threshold;
    elseif constraint == "conc"
        disp("Check Concave Constraint");
        v = diff(coef, 2) > -threshold;
    elseif constraint == "peak"
        disp("Check Peak Constraint");
        [peak, peakidx] = max(y);
        [peaksplinevalue, peaksplineidx] = max(X(peakidx, :));
        v = [diff(coef(1:peaksplineidx)) < threshold; 0; diff(coef(peaksplineidx+1:end)) > threshold];
    elseif constraint == "valley"
        disp("Check Valley Constraint");
        [valley, valleyidx] = max(-y);
        [valleysplinevalue, valleysplineidx] = max(X(valleyidx, :));
        v = [diff(coef(1:valleysplineidx)) > threshold; 0; diff(coef(valleysplineidx+1:end)) < threshold];
    elseif constraint == "none"
        v = zeros(size(coef));
    else
        disp("Constraint of type -"+constraint+"- not implemented");
    end

end
