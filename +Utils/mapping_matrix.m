function D = mapping_matrix(constraint, nr_splines)
%% 
% Create the mapping matrix for the the constraint P-spline as in Fahrmeir
% Regression, p. 436f, for the given constraint
%
% Parameters:
% -----------
% constraint : str    - Constraint type.
% nsplines   : int    - Number of used B-spline basis functions.
%
% Returns:
% --------
% D  : matrix       - Mapping matrix D
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
    constraint (1,1) string
    nr_splines (1,1) double = 10;
end
    
    if ismember(constraint, ["inc", "dec", "peak", "valley"])
        order = 1;
        d = [-1*ones(nr_splines,1), ones(nr_splines,1)];
        D = spdiags(d, 0:order, nr_splines-order, nr_splines);
    elseif ismember(constraint, ["conc", "conv", "smooth"])
        order = 2;
        d = [ones(nr_splines, 1), -2*ones(nr_splines, 1), ones(nr_splines,1)];
        D = spdiags(d, 0:order, nr_splines-order, nr_splines);
    elseif ismember(constraint, ["none"])
        D = zeros(nr_splines, nr_splines);
    else
        D = 0;
        disp("Constraint not implemented")
        return
    end

end
