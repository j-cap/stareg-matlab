function s = predict(Xpred, model, coef)
%%
% Calculate predictions for data in Xpred given the model and coef
%
%
arguments
   Xpred (:,:) double;
   model (1,1) struct;
   coef (1,1) struct;
end
    i = 1;
    fn = fieldnames(model);
    [nprediction, ~] = size(Xpred);
    B = [];

    for i=1:numel(fn)

        type_ = model.(fn{i}).type;
        nr_coef = model.(fn{i}).nr_splines; % get coefficients of the submodel
        knots = model.(fn{i}).knots;     % get knots of the submodel
        order = model.(fn{i}).spline_order; % get spline order of the submodel

        if type_.startsWith("s") % 
            dim = str2double(type_{1}(3));
            basis = zeros(nprediction, nr_coef);
            for s_idx=order.so1+1:length(knots.k1)-1 
                basis(:,s_idx-order.so1) = Bspline.basisfunction(Xpred(:,dim), knots.k1, s_idx, order.so1);
            end

        elseif type_.startsWith("t") % check the constraint for tensor-product B-splines
            dim1 = str2double(type_{1}(3));
            dim2 = str2double(type_{1}(5));
            B1 = zeros(nprediction, length(knots.k1)-1-order.so1); %basis matrix for direction 1
            B2 = zeros(nprediction, length(knots.k2)-1-order.so2); %basis matrix for direction 2
            basis = zeros(nprediction, prod(nr_coef));

            for s_idx1=order.so1+1:length(knots.k1)-1
                B1(:,s_idx1-order.so1) = Bspline.basisfunction(Xpred(:,dim1), knots.k1, s_idx1, order.so1);
            end
            for s_idx2=order.so2+1:length(knots.k2)-1
                B2(:,s_idx2-order.so2) = Bspline.basisfunction(Xpred(:,dim2), knots.k2, s_idx2, order.so2);
            end
            for data_idx=1:nprediction
                basis(data_idx,:) = kron(B2(data_idx,:), B1(data_idx,:));
            end
        else
            disp("Only B-splines (s) and tensor-product B-splines (t) available");
        end
        B = [B, basis];
        i = i + prod(nr_coef); 
    end

    % get the last coef vector out of the coef struct
    s = B*coef.("c"+string(length(fieldnames(coef))-1));

end

