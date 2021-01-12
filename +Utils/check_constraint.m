function v = check_constraint(coef, constraint, y, B)
%% 
% Check whether the coefficients in coef hold true to the constraint for
% the B-spline coefficients.
%
% Parameters:
% -----------
% coef  : array     - Array of coefficients.
% constraint : str  - Constraint type.
% y  : array        - Target data.
% B  : matrix       - B-spline basis matrix.
%
% Returns:
% --------
% v  : array      - Diagonal elements of the weighting matrix V.

arguments
   coef (:,1) double
   constraint (1,1) string = "inc"
   y (:,1) double = 0
   B (:,:) double = 0
end

    threshold = 1e-4;
    if ~ismember(constraint, ["inc", "dec", "peak","valley","conc","conv","none"])
        disp("Constraint of type -"+constraint+"- not implemented");
        return
    end
    
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
        v = [diff(coef(1:peaksplineidx)) < threshold; 0; diff(coef(peaksplineidx+1:end)) > -threshold];
    elseif constraint == "valley"
        [valley, valleyidx] = max(-y);
        [valleysplinevalue, valleysplineidx] = max(B(valleyidx, :));
        v = [diff(coef(1:valleysplineidx)) > -threshold; 0; diff(coef(valleysplineidx+1:end)) < threshold];
    elseif constraint == "none"
        v = zeros(size(coef));
    end
    v = double(v);
end
