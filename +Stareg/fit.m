function [coef, B, model] = fit(description, X,y)
%%
% Fit a structured additive regression model using B-splines and
% tensor-product B-splines to the data given in (X, y) following the
% description
% 
% Parameters:
% -----------
% description : struct with fields    - Describes the model structure, e.g. 
%                       create the cell array description of the model

%                       d = {["s(1)", "inc", 100, 3000, "e"]; 
%                            ["t(1,2)", "inc,none", "12,20", "2000,2000", "e,q"]};
%                       describing a model using a P-spline with increasing constraint and 100 
%                       basis functions for dimension 1 and a tensor-product P-spline with
%                       increasing constraint for dimension 1 and no constraints for dimension 2
%                       using 12 and 10 basis functions for the respective dimension.
%
% X :  array        - input data, shape (n_samples, n_dim)
% y : array         - target data, shape (n_samples, )
%    
% Returns:
% coef : array              - Estimated coefficients
% basis_matrix : matrix     - Basis matrix to get evaluations
% model : struct            - Model struct
%   

arguments
    description (:,1) cell;
    X (:,:) double;
    y (:,1) double;
end

    model = Stareg.create_model_from_description(description, X, y);
    [B, S, K, W, coef_pls] = Stareg.create_model_matrices(model);
    coef = struct("c0", coef_pls);

    iterIdx = 1;
    BtB = B'*B;
    Bty = B'*y;

    W_compare = W;
    W_old = zeros(size(W));

    while ~all(W_old == W_compare)
        W_old = W_compare;
        coef_pls = (BtB + S + K) \ Bty;
        coef.("c"+string(iterIdx)) = coef_pls;
        
        [w, W_compare] = Utils.check_constraint_full(coef_pls, model, B, y);
        K = Utils.create_constraint_matrix(model, w);
        
        fprintf("Iteration %d \n", iterIdx);
        fprintf("MSE = %4.4f \n", Utils.mse(y, B*coef_pls));
        fprintf("---------------------------------------\n");

        iterIdx = iterIdx + 1;
        if iterIdx > 25
            disp("Stop the count")
            return
        end
    end
    
end

