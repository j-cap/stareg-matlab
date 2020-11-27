function D = mapping_matrix(constraint, nr_splines)
%% Create the mapping matrix for the given constraint.
%
% Parameters:
% -----------
% constraint : str    - Constraint type.
% nsplines   : int    - Number of splines.
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

e = ones(nr_splines, 1);
if constraint == "inc" | constraint == "dec" | constraint == "peak" | constraint == "valley"
    diagonals = [-e, e];
    degree = 1;
elseif constraint == "conv" | constraint == "conc" | constraint == "smooth"
    diagonals = [e, -2*e, e];
    degree = 2;
elseif constraint == "none"
    diagonals = zeros(nr_splines, 1);
    degree = 0;
else
    disp("Constraint of type -"+constraint+"- not implemented");
    return
end

D = spdiags(diagonals, 0:degree, nr_splines-degree, nr_splines);

end
