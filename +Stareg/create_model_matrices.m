function [basis_matrix, smoothness_matrix, constr_matrix, weight_matrix, coef] = create_model_matrices(model)
%%
% create basis matrix, smoothness matrix and constraint matrix for the
% given model
%
% Parameters:
% -----------
% model : struct with fields    - Created by the function
%                                 Stareg.create_model_from_description
%
% Returns:
% --------
% basis_matrix : matrix
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

