function s = storage_efficient_prediction(Xpred, reduced_model, coef)
%%
% Calculate predictions for data in Xpred given the model and coef.
%
% Inputs:
% -------
% Xpred : matrix           - Input data to evaluate the model on
% reduced_model : struct   - Minimal Model struct.
% coef : array             - Coefficients to calculate the prediction for
%
% Outputs:
% --------
% s : array            - Predictions for input data Xpred.
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
   Xpred (:,:) double;
   reduced_model (1,1) struct;
   coef (:,1) double;
end

    fn = fieldnames(reduced_model);
    [n_pred, ~] = size(Xpred);
    B = [];
 
    for submodel=1:numel(fn) % iterate over all submodels

        type_ = reduced_model.(fn{submodel}).type;
        nr_coef = reduced_model.(fn{submodel}).nr_splines; % get coefficients of the submodel
        knots = reduced_model.(fn{submodel}).knots;     % get knots of the submodel
        order = reduced_model.(fn{submodel}).spline_order; % get spline order of the submodel

        if type_.startsWith("s") % 
            dim = str2double(type_{1}(3));
            basis = zeros(n_pred, nr_coef);
            for spline_idx=order.so1+1:length(knots.k1)-1 
                basis(:,spline_idx-order.so1) = Bspline.basisfunction(Xpred(:,dim), knots.k1, spline_idx, order.so1);
            end

        elseif type_.startsWith("t") % check the constraint for tensor-product B-splines
            dim1 = str2double(type_{1}(3));
            dim2 = str2double(type_{1}(5));
            B1 = zeros(n_pred, length(knots.k1)-1-order.so1); %basis matrix for direction 1
            B2 = zeros(n_pred, length(knots.k2)-1-order.so2); %basis matrix for direction 2
            basis = zeros(n_pred, prod(nr_coef));

            for spline_idx1=order.so1+1:length(knots.k1)-1
                B1(:,spline_idx1-order.so1) = Bspline.basisfunction(Xpred(:,dim1), knots.k1, spline_idx1, order.so1);
            end
            for spline_idx2=order.so2+1:length(knots.k2)-1
                B2(:,spline_idx2-order.so2) = Bspline.basisfunction(Xpred(:,dim2), knots.k2, spline_idx2, order.so2);
            end
            for data_idx=1:n_pred
                basis(data_idx,:) = kron(B2(data_idx,:), B1(data_idx,:));
            end
        else
            disp("Only B-splines (s) and tensor-product B-splines (t) available");
        end
        B = [B, basis];
    end
    
    s = B*coef;

end

