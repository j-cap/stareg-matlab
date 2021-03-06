function [coef_, B, model, reduced_model, coef_list] = fit(description, X,y)
%%
% Fit a structured additive regression model using B-splines and
% tensor-product B-splines to the data given in (X, y) following the
% description
% 
% Inputs:
% -----------
% description : struct with fields    
%                   - Describes the model structure, e.g. 
%                       create the cell array description of the model
%                       d = {["s(1)", 100, "inc", 3000, "e"]; 
%                            ["t(1,2)", "12,20", "inc,none", "2000,2000", "e,q"]};
%                       describing a model using a P-spline with increasing constraint and 100 
%                       basis functions for dimension 1 and a tensor-product P-spline with
%                       increasing constraint for dimension 1 and no constraints for dimension 2
%                       using 12 and 10 basis functions for the respective dimension.
%
% X : array         - input data, shape (n_samples, n_dim)
% y : array         - target data, shape (n_samples, )
%    
% Outputs:
% --------
% coef_ : array             - Estimated coefficients.
% B : matrix                - Basis matrix to evaluation the model. 
% model : struct            - Model struct.
% reduced_model : struct    - Model struct for fast_prediction.
% coef_list : struct        - List of the coeffients over all iterations.
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
    description (:,1) cell;
    X (:,:) double;
    y (:,1) double;
end

    model = Stareg.create_model_from_description(description, X, y); % create model from description
    [B, S, K, V, coef_pls] = Stareg.create_model_matrices(model); % create model matrices
    coef_list = struct("c0", coef_pls); %  initialized coef struct for the iteration

    iterIdx = 1;
    BtB = B'*B; % pre-calculate basis matrix product
    Bty = B'*y; % pre-calculate basis matrix and target data product

    V_compare = V;
    V_old = zeros(size(V)); % create initial weight matrix for comparison

    while ~all(V_old == V_compare) % check condition for termination
        V_old = V_compare;
        coef_pls = (BtB + S + K) \ Bty; % calculate new coefficients
        coef_list.("c"+string(iterIdx)) = coef_pls; % add new coefficients to coef struct
        
        [v, V_compare] = Utils.check_constraint_full(coef_pls, model, B, y); % check the new coef against all constraints
        K = Utils.create_constraint_matrix(model, v); % create new constraint matrix K
        
        % print some metrics
        fprintf("Iteration %d \n", iterIdx);
        fprintf("MSE = %4.4f \n", Utils.mse(y, B*coef_pls));
        fprintf("---------------------------------------\n");

        iterIdx = iterIdx + 1;
        if iterIdx > 25
            disp("Stop the count")
            break
        end
    end
    % create reduced model for fast prediciton
    coef_ = coef_pls;
    [reduced_model] = Stareg.create_reduced_model(model);
end

