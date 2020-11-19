function D = mapping_matrix_tps(constraint, nsplines)
%%
% Create the mapping matrix for the given constraint for tensorproduct splines.
% In contrast to the 1D case and B-splines, we now return the matrix 
%               D = sqrt(D_c'*D_c)
% such that K = D'VD can be computed in PIRLS
% 
% Parameters:
% -----------
% constraint : str  - constraint type
% nsplines : vector    - Nr. of splines for dim 1 and 2.
%
% Returns:
% --------
% D  : matrix       - Mapping matrix D
%
%%
if constraint == "inc" | constraint == "dec" | constraint == "peak" | constraint == "valley"
    D = sqrtm(full(smoothness_matrix_TPS(nsplines(1), nsplines(2), 1)));
elseif constraint == "conv" | constraint == "conc"
    D = sqrtm(full(smoothness_matrix_TPS(nsplines(1), nsplines(2), 2)));
elseif constraint == "none"
    D = eye(prod(ns));
else
    disp("Constraint of type -"+constraint+"- not implemented");
end
end
