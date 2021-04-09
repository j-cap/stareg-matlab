function [v, vc] = check_constraint_full(coef, model, basis, y)
%%
% Checks the respective parts of the coef vector against the respective
% constraint 
%
% Parameters:
% -----------
% coef : array         - Vector of coefficients.
% model : struct       - Model structure.
% basis : matrix 
% y : array            - Target variable.
%
% Returns:
% --------
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
    coef (:,1) double;
    model (1,1) struct;
    basis (:,:) double;
    y (:,1) double;
end


    vc = [];
    v = struct();
    coef_idx = 1;
    fn = fieldnames(model);
    
    for i=1:numel(fn)

        nr_coef = model.(fn{i}).nr_splines; % get coefficients of the submodel
        constr = model.(fn{i}).constraint; %  get constraint of the submodel

        if model.(fn{i}).type.startsWith("s") % check the constraint for B-splines
            test_coef = coef(coef_idx:coef_idx+nr_coef-1);
            vv = Utils.check_constraint(test_coef, constr, y, basis(:,[i:i+nr_coef-1]));
            vc = [vc; vv];
            v.(fn{i}) = struct("v", vv);
        elseif model.(fn{i}).type.startsWith("t") % check the constraint for tensor-product B-splines
            test_coef = coef(coef_idx:coef_idx+prod(nr_coef)-1);
            v1 = Utils.check_constraint_dim1(test_coef, constr(1), nr_coef);
            v2 = Utils.check_constraint_dim2(test_coef, constr(2), nr_coef);
            vc = [vc; v1; v2];
            v.(fn{i}) = struct("v1", v1, "v2", v2);
        else
            disp("Only B-splines (s) and tensor-product B-splines (t) available");
        end
        coef_idx = coef_idx + prod(nr_coef); 

    end
end

