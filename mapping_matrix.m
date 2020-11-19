function D = mapping_matrix(constraint, nsplines, knots)
%% Create the mapping matrix for the given constraint.
%
% Parameters:
% -----------
% constraint : str    - Constraint type.
% nsplines   : int    - Number of splines.
% knots      : array  - Either knots or 0. If 0, we use the FD matrix.
%
% Returns:
% --------
% D  : matrix       - Mapping matrix D
%
%%
if constraint == "inc" | constraint == "dec" | constraint == "peak" | constraint == "valley"
    D = smoothness_matrix(nsplines, 1, knots);
elseif constraint == "conv" | constraint == "conc"
    D = smoothness_matrix(nsplines, 2, knots);
elseif constraint == "none"
    D = eye(nsplines);
else
    disp("Constraint of type -"+constraint+"- not implemented");
end

end
