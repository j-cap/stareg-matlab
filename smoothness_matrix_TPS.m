function S = smoothness_matrix_TPS(nparam1,nparam2,degree)
%%
%
% Parameters:
% ------------
%   nparam1 : int  - Number of splines for dimension 1
%   nparam2 : int  - Number of splines for dimension 2
%   degree  : int  - Order of the smoothness matrix (1 or 2)
%
% Returns:
% ------------
%   S     : matrix - Smoothness penalty matrix
%%
    S1 = smoothness_matrix(nparam1, degree, 0);
    S2 = smoothness_matrix(nparam2, degree, 0);
    K1 = S1'*S1;
    K2 = S2'*S2;
    S = kron(eye(nparam2), K1) + kron(K2, eye(nparam1));
end

