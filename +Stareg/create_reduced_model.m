function [reduced_model] = create_reduced_model(model)
%%
% Creates the reduced_model for fast_prediction based on the given model struct.
%

% Inputs:
% -----------
% model : struct            - Complete model description.
%
% Outputs:
% --------
% reduced_model : struct    - Reduced model for fast_prediction.
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
end
    
    reduced_model = struct;
    fn = fieldnames(model);
       
    for submodel=1:numel(fn)
        reduced_model.(string(fn(submodel))) = struct();
    end

    for k=1:numel(fn) % iterate over all submodels
        submodel = model.(fn{k}); % get the submodel
        reduced_model.(fn{k}).type = submodel.type;
        reduced_model.(fn{k}).nr_splines = submodel.nr_splines;
        reduced_model.(fn{k}).knots = submodel.knots;
        reduced_model.(fn{k}).spline_order = submodel.spline_order;
    end

end

