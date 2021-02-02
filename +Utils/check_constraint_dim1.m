function v1 = check_constraint_dim1(coef, constraint, nr_splines)
%%
% Compute the diagonal elements of the weighting matrix for SC-TP-P-splines 
% given the constraint for direction 1.
%    
% According to the scheme given in the Master Thesis !!
%
% Parameters:
% -----------
% coef  : array      - Coefficient vector to test against constraint.
% constraint : str   - Specifies the constraint.
% nr_splines : list  - Specifies the number of splines in each dimension
%
% Returns
% -------
% v1  : array         - Diagonal elements of the weighting matrix V.
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
    constraint (1,1) string = "inc";
    nr_splines (1,2) double = [12,8];
end
    switch constraint
        case "inc"
            diff = 1;
        case "dec"
           diff = 1;
        otherwise
            diff = 0;
    end

    v1 = zeros((nr_splines(1)-diff)*nr_splines(2), 1);
    for i=1:nr_splines(2)
        for j=1:nr_splines(1)-1
            v1(j+(i-1)*(nr_splines(1)-1)) = coef(j+(i-1)*nr_splines(1) + 1) - coef(j+(i-1)*nr_splines(1));
        end
    end

    switch constraint
        case "inc"
            v1 = v1 < 0;
        case "dec"
            v1 = v1 > 0;
        case "none"
            v1 = zeros(size(v1));
        otherwise
            disp("Constraint not implemented!");
            return
    end

end

