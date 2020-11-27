function v = check_constraint(coef, constraint, B, y)
%% 
% Check whether the coefficients in coef hold the constraint
%
% Parameters:
% -----------
% coef  : array     - Array of coefficients.
% constraint : str  - Constraint type.
% B  : matrix       - B-spline basis matrix.
% y  : array        - Output data.
%
% Returns:
% --------
% v  : array      - Diagonal elements of the weighting matrix V.

arguments
   coef (:,1) double
   constraint (1,1) string = "inc"
   B (:,:) double = 0
   y (:,1) double = 0
end

    threshold = 0;
    if constraint == "inc"
        v = diff(coef) < threshold;
    elseif constraint == "dec"
        v = diff(coef) > -threshold;
    elseif constraint == "conv"
        v = diff(coef, 2) < threshold;
    elseif constraint == "conc"
        v = diff(coef, 2) > -threshold;
    elseif constraint == "peak"
        [peak, peakidx] = max(y);
        [peaksplinevalue, peaksplineidx] = max(B(peakidx, :));
        v = [diff(coef(1:peaksplineidx)) < threshold; 0; diff(coef(peaksplineidx+1:end)) > threshold];
    elseif constraint == "valley"
        [valley, valleyidx] = max(-y);
        [valleysplinevalue, valleysplineidx] = max(B(valleyidx, :));
        v = [diff(coef(1:valleysplineidx)) > threshold; 0; diff(coef(valleysplineidx+1:end)) < threshold];
    elseif constraint == "none"
        v = zeros(size(coef));
    else
        disp("Constraint of type -"+constraint+"- not implemented");
    end

end
