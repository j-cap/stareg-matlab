function [constr_matrix] = create_constraint_matrix(model, weights)
%%
% Creates the constraint matrix using the new weights in weights
%
% Parameters:
% -----------
% model : struct with fields
% weights : struct with fields
% 
% Returns:
% --------
% K : matrix
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
    model (1,1) struct;
    weights (1,1) struct;
end
    fn = fieldnames(model);
    constr_matrix = [];
    for i=1:numel(fn)
        if model.(fn{i}).type.startsWith("s") % create smoothness matrix and constraint matrix for B-splines
            constr_matrix = blkdiag(constr_matrix, model.(fn{i}).lam_c * model.(fn{i}).Dc' * diag(weights.(fn{i}).v) * model.(fn{i}).Dc);

        elseif model.(fn{i}).type.startsWith("t") % create smoothness matrix and constraint matrix for tensor-product B-splines

            Dc1 = model.(fn{i}).lam_c(1) * (model.(fn{i}).Dc.Dc1' * diag(weights.(fn{i}).v1) * model.(fn{i}).Dc.Dc1);
            Dc2 = model.(fn{i}).lam_c(2) * (model.(fn{i}).Dc.Dc2' * diag(weights.(fn{i}).v2) * model.(fn{i}).Dc.Dc2);
            Dc = Dc1 + Dc2;
            constr_matrix = blkdiag(constr_matrix, Dc);
                        
        end
    end

end

