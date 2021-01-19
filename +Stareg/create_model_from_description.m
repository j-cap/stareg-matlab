function [model] = create_model_from_description(description,X, y)
%%
% Create the model struct with fields from the given description  
%    
% Parameters:
% -----------
% descr : tuple of tuples    - Describes the model structure, e.g. 
%                              %% create the cell array description of the model
%                                 d = {["s(1)", "inc", 100, 3000, "e"]; 
%                                      ["t(1,2)", "none,none", "12,20", "2000,2000", "e,q"]};
%                                 describing a model using a P-spline with increasing constraint and 100 
%                                 basis functions for dimension 1 and a tensor-product P-spline without 
%                                 constraints using 12 and 10 basis functions for the respective dimension.
% X : array                  - np.array of the input data, shape (n_samples, n_dim)
% y : array                  - np.array of the target data, shape (n_samples, )
%    
% Returns:
% --------
% d : dict      - Returns a dictionary with the following key-value pairs:
%                        basis=B, 
%                        smoothness=S, 
%                        constraint=K, 
%                        opt_lambdas=optimal_lambdas, 
%                        coef_=coef_pls, 
%                        weights=weights

arguments
   description (:,1) cell;
   X (:,:) double;
   y (:,1) double;
end

%%      
    model_base = struct;
    model = struct;
    for i=1:numel(description)
        model_base.("f"+string(i)) = description{i}; 
    end

    fn = fieldnames(model_base);
    for k=1:numel(fn)
       s = model_base.(fn{k});
       type_ = s(1);
       constr = s(2);
       nr_splines = s(3);
       lam_c = s(4);
       knot_type = s(5);
       if s(1).startsWith("s")
           dim = str2double(type_{1}(3));
           nr_splines = str2double(nr_splines);
           lam_c = str2double(lam_c);
           [B,~] = Bspline.basismatrix(X(:,dim), nr_splines, 3, knot_type);
           Ds = Utils.mapping_matrix("smooth", nr_splines);
           Dc = Utils.mapping_matrix(constr, nr_splines);
           best_lam = Bspline.calc_GCV(X(:,dim), y, 100, 0, nr_splines, 3, knot_type);
           [coef_pls, B_, k_] = Bspline.fit_Pspline(X(:,dim),y,best_lam, nr_splines, 3, knot_type);
           v = Utils.check_constraint(coef_pls, constr, y, B);
           model.(fn{k}) = struct("type", type_, "B",B,"Ds",Ds,"Dc",Dc,"best_lam",best_lam,"lam_c", lam_c, "coef_pls",coef_pls, "v", v, ...
               "knot_type", knot_type, "spline_order", 3, "nr_splines", nr_splines, "constraint", constr);
       elseif s(1).startsWith("t")
           dim1 = str2double(type_{1}(3));
           dim2 = str2double(type_{1}(5));
           constr = constr(1).split(",");
           nr_splines = str2double(nr_splines.split(","));
           lam_c = str2double(lam_c.split(","));
           knot_type = knot_type.split(",");
           [T,~] = Bspline.tensorproduct_basismatrix(X(:,[dim1,dim2]), nr_splines, [3,3], knot_type);
           Ds1 = Utils.mapping_matrix_tp("smooth", nr_splines, 1);
           Ds2 = Utils.mapping_matrix_tp("smooth", nr_splines, 2);
           Dc1 = Utils.mapping_matrix_tp(constr(1), nr_splines, 1);
           Dc2 = Utils.mapping_matrix_tp(constr(2), nr_splines, 2);
           best_lam = Bspline.calc_GCV_2d(X(:,[dim1,dim2]), y, 100, 0, nr_splines, [3,3], knot_type);
           [coef_pls, T_, k_] = Bspline.fit_Pspline(X(:,[dim1,dim2]), y, best_lam, nr_splines, [3,3], knot_type);
           v1 = Utils.check_constraint_dim1(coef_pls, constr(1), nr_splines);
           v2 = Utils.check_constraint_dim2(coef_pls, constr(2), nr_splines);
           model.(fn{k}) = struct("type", type_, "B",T,"Ds",struct("Ds1", Ds1, "Ds2", Ds2),"Dc", struct("Dc1", Dc1, "Dc2", Dc2), ...
               "best_lam",best_lam,"lam_c", lam_c, "coef_pls",coef_pls, "v", struct("v1", v1, "v2", v2), "knot_type", knot_type, ...
               "spline_order", struct("so1", 3, "so2", 3), "nr_splines", nr_splines, "constraint", constr);
       end
    end
    
    
end

