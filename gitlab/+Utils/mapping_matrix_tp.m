function Dc = mapping_matrix_tp(constraint, nr_splines, dim)
%%
% Create the mapping matrix for the constraint tensorproduct P-splines as 
% in Fahrmeir, Regression, p.508 (8.27) for the constraint.
%
% Parameters:
% -----------
% constraint : str    - Constraint type.
% nsplines   : array  - Number of splines for both dimensions.
% dim : int           - indicator for the dimension of the constraint, 
%                       1 for dimension 1, 2 for dimension 2, e.g. 
%                       (10, "inc", 1) == 10 basis functions with
%                       increasing constraint in dimension 1 for the
%                       two-dimesional data X = [x_1, x_2].
% Returns:
% --------
% Dc : matrix       - Mapping matrix for the constraint and dimension.
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
    nr_splines (1,2) double = [10,10];
    dim (1,1) double = 1;
end
    if ismember(constraint, ["inc", "dec", "peak", "valley"]) % for finite difference order 1 constraints
        order = 1;
        d = [-1*ones(nr_splines(dim),1), ones(nr_splines(dim),1)];
        D = spdiags(d, 0:order, nr_splines(dim)-order, nr_splines(dim));
    elseif ismember(constraint, ["conc", "conv", "smooth"]) % for finited difference order 2 constraints
        order = 2;
        d = [ones(nr_splines(dim), 1), -2*ones(nr_splines(dim), 1), ones(nr_splines(dim),1)];
        D = spdiags(d, 0:order, nr_splines(dim)-order, nr_splines(dim));
    elseif constraint == "none"
        Dc = zeros(prod(nr_splines));
        return
    end
    
    if dim == 1
        Dc = kron(eye(nr_splines(dim+1)), D);
    elseif dim == 2
        Dc = kron(D, eye(nr_splines(dim-1)));
    end 
end
