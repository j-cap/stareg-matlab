function v2 = check_constraint_dim2(coef, constraint, nr_splines)
%%
% Compute the diagonal elements of the weighting matrix for SC-TP-P-splines 
% given the constraint for direction 2.
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
% v2  : array         - Diagonal elements of the weighting matrix V.
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
arguments
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

    v2 = zeros(nr_splines(1)*(nr_splines(2)-diff), 1);
    for i=1:nr_splines(2)-1
        for j=1:nr_splines(1)
            v2(j+(i-1)*(nr_splines(1))) = coef(j+i*nr_splines(1)) - coef(j+i*nr_splines(1)-nr_splines(1));
        end
    end

    switch constraint
        case "inc"
            v2 = v2 < 0;
        case "dec"
            v2 = v2 > 0;
        case "none"
            v2 = zeros(size(v2));
        otherwise
            disp("Constraint not implemented!");
            return
    end

end

