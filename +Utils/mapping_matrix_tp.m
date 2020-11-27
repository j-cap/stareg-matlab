function [Id2D1, D2Id1] = mapping_matrix_tp(constraint, nr_splines)
%% Create the mapping matrix for the given constraint.
%
% Parameters:
% -----------
% constraint : str    - Constraint type.
% nsplines   : array  - Number of splines for both dimensions.
%
% Returns:
% --------
% Id2D1 : matrix       - Mapping matrix for dimension 1.
% D2Id1 : matrix       - Mapping matrix for dimension 2.
%%
arguments
    constraint (1,1) string
    nr_splines (1,2) double = [10,10]
end

e1 = ones(nr_splines(1), 1);
e2 = ones(nr_splines(2), 1);
if constraint == "inc" | constraint == "dec" | constraint == "peak" | constraint == "valley"
    d1 = [-e1, e1];
    d2 = [-e2, e2];
    degree = 1;
elseif constraint == "conv" | constraint == "conc" | constraint == "smooth"
    d1 = [e1, -2*e1, e1];
    d2 = [e2, -2*e2, e2];
    degree = 2;
elseif constraint == "none"
    d1 = zeros(nr_splines(1), 1);
    d2 = zeros(nr_splines(2), 1);
    degree = 0;
else
    disp("Constraint of type -"+constraint+"- not implemented");
    return
end

D1 = spdiags(d1, 0:degree, nr_splines(1)-degree, nr_splines(1));
D2 = spdiags(d2, 0:degree, nr_splines(2)-degree, nr_splines(2));

Id2D1 = kron(eye(nr_splines(2)), D1);
D2Id1 = kron(D2, eye(nr_splines(1)));

end
