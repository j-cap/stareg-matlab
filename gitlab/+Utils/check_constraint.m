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
%
% Dependencies:
%    Matlab release: 2020b
%
% This function is part of: stareg-matlab
%
% Author:  Jakob Weber
% email:   jakob.weber@ait.ac.at
% Company: Austrian Institute of Technology GmbH
%          Complex Dynamical Systems
%          Center for Vision, Automation & Control
%          http://www.ait.ac.at
%
% Version: 1.0.1 - 2021-02-02

% Change log:
% x.y.z - 2021-02-02 - author:
% - added important feature, s. issue #34
% - fixed bug #2
%%
arguments % default values
   coef (:,1) double
   constraint (1,1) string = "inc"
   y (:,1) double = 0
   B (:,:) double = 0
end

    threshold = 1e-6;
    if ~ismember(constraint, ["inc", "dec", "peak","valley","conc","conv","none"])
        disp("Constraint of type -"+constraint+"- not implemented");
        return
    end
    
    switch constraint
        case "inc"
            v = diff(coef) < threshold;
        case "dec"
            v = diff(coef) > -threshold;
        case "conv"
            v = diff(coef, 2) < threshold;
        case "conc"
            v = diff(coef, 2) > -threshold;
        case "peak"
            [peak, peakidx] = max(y);
            [peaksplinevalue, peaksplineidx] = max(B(peakidx, :));
            v = [diff(coef(1:peaksplineidx)) < threshold; 0; diff(coef(peaksplineidx+1:end)) > -threshold];
        case "valley"
            [valley, valleyidx] = max(-y);
            [valleysplinevalue, valleysplineidx] = max(B(valleyidx, :));
            v = [diff(coef(1:valleysplineidx)) > -threshold; 0; diff(coef(valleysplineidx+1:end)) < threshold];
        case "none"
            v = zeros(size(coef));
        otherwise
            disp("Constraint not implemented!");
            return
    end
    v = double(v);
end
