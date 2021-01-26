function [basis_matrix, smoothness_matrix, constr_matrix, weight_matrix, coef] = create_model_matrices(model)
%%
% Creates combined basis matrix, smoothness matrix and constraint matrix for the
% given model, i.e.
%
% for the model given by the description 
%    description = {
%        ["s(1)", "inc", 100, 3000, "e"]; 
%        ["s(2)", "peak", 10, 300, "e"];
%        ["t(1,2)", "inc,none", "12,20", "2000,2000", "e,q"]
%        };
% we generate 
%    - the basis_matrix B = [B_1, B_2, B_3] as horizontal concatenation 
%    - the smoothness_matrix S = blkdiag(S_1, S_2, S_3), with the
%       individal smoothness matrices, e.g. S_1 = \lamda_1 * D_2'*D_2, etc.
%    - the constraint_matrix K = blkdiag(K_1, K_2, K_3), with the
%       individual constraint matrice, e.g. K_1 = \lambda_c1 * D_c'*D_c, etc.
%    - the weight_matrix W = diag(w_1, w_2, w_31, w_32) as diagonal matrix
%       of the individual weights w_1, w_2, w_31, w_w32 (w_31 means weights
%       for the 3rd function aka. t(1,2) and constraint 1)
%    - the coefficient vector coef = [coef_1, coef_2, coef_3] as vertical
%       concatenation of the individual coefficient vectors
%
%
% Parameters:
% -----------
% model : struct with fields    - Created by the function
%                                 Stareg.create_model_from_description
%
% Returns:
% --------
% basis_matrix : matrix         - basis_matrix 
% smoothness_matrix : matrix
% constr_matrix : matrix    
% weight_matrix : matrix           - IMPORANT: only for comparison reason,
%                                    does not include the original weights
%                                    for tensor-product B-splines.
% coef : array 
%%
arguments
    model (1,1) struct;
end
    fn = fieldnames(model);
    basis_matrix = [];
    constr_matrix = [];
    smoothness_matrix = [];
    weight_matrix = [];
    coef = [];
    for i=1:numel(fn)

        basis_matrix = [basis_matrix, model.(fn{i}).B]; % create basis
        if model.(fn{i}).type.startsWith("s") % create smoothness matrix and constraint matrix for B-splines
            smoothness_matrix = blkdiag(smoothness_matrix, model.(fn{i}).best_lam * model.(fn{i}).Ds' * model.(fn{i}).Ds);
            constr_matrix = blkdiag(constr_matrix, model.(fn{i}).lam_c * model.(fn{i}).Dc' * diag(model.(fn{i}).v) * model.(fn{i}).Dc);
            weight_matrix = [weight_matrix; (model.(fn{i}).v)];
        elseif model.(fn{i}).type.startsWith("t") % create smoothness matrix and constraint matrix for tensor-product B-splines
            Ds = model.(fn{i}).best_lam * ( model.(fn{i}).Ds.Ds1' * model.(fn{i}).Ds.Ds1 + model.(fn{i}).Ds.Ds2' * model.(fn{i}).Ds.Ds2);
            smoothness_matrix = blkdiag(smoothness_matrix, Ds);

            Dc1 = model.(fn{i}).lam_c(1) * (model.(fn{i}).Dc.Dc1' * diag(model.(fn{i}).v.v1) * model.(fn{i}).Dc.Dc1);
            Dc2 = model.(fn{i}).lam_c(2) * (model.(fn{i}).Dc.Dc2' * diag(model.(fn{i}).v.v2) * model.(fn{i}).Dc.Dc2);
            Dc = Dc1 + Dc2;
            constr_matrix = blkdiag(constr_matrix, Dc);
            
            weight_matrix = [weight_matrix; model.(fn{i}).v.v1; model.(fn{i}).v.v2];
            
        end
        coef = [coef; model.(fn{i}).coef_pls];
    end

end

