function D = mapping_matrix(constraint, nr_splines)
%% Create the mapping matrix for the the constraint P-spline as in Fahrmeir
% Regression, p. 436f, for the given constraint
%
% Parameters:
% -----------
% constraint : str    - Constraint type.
% nsplines   : int    - Number of used B-spline basis functions.
%
% Returns:
% --------
% D  : matrix       - Mapping matrix D
%
%%
arguments
    constraint (1,1) string
    nr_splines (1,1) double = 10;
end

    if ismember(constraint, ["inc", "dec", "peak", "valley"])
        order = 1;
    else
        order = 2;
    end
    
    if order == 1
        d = [-1*ones(nr_splines,1), ones(nr_splines,1)];
        D = spdiags(d, 0:order, nr_splines-order, nr_splines);
    elseif order == 2
        d = [ones(nr_splines, 1), -2*ones(nr_splines, 1), ones(nr_splines,1)];
        D = spdiags(d, 0:order, nr_splines-order, nr_splines);
    end

    if constraint == "none"
        D = zeros(nr_splines, nr_splines);
    end

end
