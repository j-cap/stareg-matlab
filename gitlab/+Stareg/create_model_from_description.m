function [model] = create_model_from_description(description,X, y)
%%
% Create the model struct with fields from the given description.
%
% Note: the description needs to be a cell array for which each cell needs
%       to follow the following pattern
%       {[type of spline, number of splines, constraint, constraint parameter, knot_type]}
%    
% Note: for tensor-product B-splines, the values in the cell need to be strings.
%
% Inputs:
% -----------
% descr : cell array  - Describes the model structure, e.g. 
%                       create the cell array description of the model
%                       d = {["t(1,2)", "12,20", "none,none", "2000,2000", "e,q"]
%                            ["s(1)", 100, "inc", 3000, "e"] };
%                       describing a model using a P-spline with increasing constraint and 100 
%                       basis functions for dimension 1 and a tensor-product P-spline without 
%                       constraints using 12 and 10 basis functions for the respective dimension.
% X : array          - Input data of shape (n_samples, n_dim)
% y : array          - Target data of shape (n_samples, )
%    
% Outputs:
% --------
% model : struct  - Model description with all information. 
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

%%      
    model_base = struct;
    model = struct;
    for i=1:numel(description)
        model_base.("f"+string(i)) = description{i}; 
    end

    fn = fieldnames(model_base);
    for k=1:numel(fn) % iterate over all submodels
       submodel = model_base.(fn{k}); % get the submodel
       type_ = submodel(1);
       nr_splines = submodel(2);
       constr = submodel(3);
       lam_c = submodel(4);
       knot_type = submodel(5);
       if submodel(1).startsWith("s") % create matrices for type is B-spline
           dim = str2double(type_{1}(3));
           nr_splines = str2double(nr_splines);
           lam_c = str2double(lam_c);
           [B,knots] = Bspline.basismatrix(X(:,dim), nr_splines, 3, knot_type); % create basismatrix
           Ds = Utils.mapping_matrix("smooth", nr_splines); % create mapping matrix for smoothness penalty
           Dc = Utils.mapping_matrix(constr, nr_splines); % create mapping matrix for constraint penalty
           best_lam = Bspline.calc_GCV(X(:,dim), y, 100, 0, nr_splines, 3, knot_type); % find optimal lambda using GCV
           [coef_pls, B_, k_] = Bspline.fit_Pspline(X(:,dim),y,best_lam, nr_splines, 3, knot_type); % calculate Penalized Least-Squares fit using optimal lambda
           v = Utils.check_constraint(coef_pls, constr, y, B); % get initial weight vector for the constraint
           model.(fn{k}) = struct("type", type_, "B",B,"Ds",Ds,"Dc",Dc,"best_lam",best_lam,"lam_c", lam_c, "coef_pls",coef_pls, "v", v, ...
               "knot_type", knot_type, "spline_order", struct("so1", 3), "nr_splines", nr_splines, "constraint", constr, "knots", knots);
       elseif submodel(1).startsWith("t") % create matrices for type is tensor-product B-spline
           dim1 = str2double(type_{1}(3)); 
           dim2 = str2double(type_{1}(5));
           constr = constr(1).split(",");
           nr_splines = str2double(nr_splines.split(","));
           lam_c = str2double(lam_c.split(","));
           knot_type = knot_type.split(",");
           [T,knots] = Bspline.tensorproduct_basismatrix(X(:,[dim1,dim2]), nr_splines, [3,3], knot_type); % create basismatrix
           Ds1 = Utils.mapping_matrix_tp("smooth", nr_splines, 1); % create mapping matrix for smoothness penalty and dimension 1
           Ds2 = Utils.mapping_matrix_tp("smooth", nr_splines, 2); % create mapping matrix for smoothness penalty and dimension 2
           Dc1 = Utils.mapping_matrix_tp(constr(1), nr_splines, 1); % create mapping matrix for constraint penalty and dimension 1
           Dc2 = Utils.mapping_matrix_tp(constr(2), nr_splines, 2); % create mapping matrix for constraint penalty and dimension 2
           best_lam = Bspline.calc_GCV_2d(X(:,[dim1,dim2]), y, 100, 0, nr_splines, [3,3], knot_type); % find optimal lambda using GCV
           [coef_pls, T_, k_] = Bspline.fit_Pspline(X(:,[dim1,dim2]), y, best_lam, nr_splines, [3,3], knot_type); % calculate Penalized Least-Squares fit using optimal lambda
           v1 = Utils.check_constraint_dim1(coef_pls, constr(1), nr_splines); % get initial weight vector for the constraint in dimension 1
           v2 = Utils.check_constraint_dim2(coef_pls, constr(2), nr_splines); % get initial weight vector for the constraint in dimension 2
           model.(fn{k}) = struct("type", type_, "B",T,"Ds",struct("Ds1", Ds1, "Ds2", Ds2),"Dc", struct("Dc1", Dc1, "Dc2", Dc2), ...
               "best_lam",best_lam,"lam_c", lam_c, "coef_pls",coef_pls, "v", struct("v1", v1, "v2", v2), "knot_type", knot_type, ...
               "spline_order", struct("so1", 3, "so2", 3), "nr_splines", nr_splines, "constraint", constr,"knots", knots);
       end
    end
    
    
end

