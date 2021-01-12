function v1 = check_constraint_dim1(coef, constraint, nr_splines)
%%
% Compute the diagonal elements of the weighting matrix for SC-TP-P-splines 
% given the constraint for direction 1.
%    
% According to the scheme given in the Master Thesis !!
%
% Parameters:
% -----------
% coef  : array      - Coefficient vector to test against constraint.
% constraint : str   - Specifies the constraint.
% nr_splines : list  - Specifies the number of splines in each dimension
%
% Returns
% -------
% v1  : array         - Diagonal elements of the weighting matrix V.
%%
arguments
    coef (:,1) double
    constraint (1,1) string = "inc";
    nr_splines (1,2) double = [12,8];
end
    if ismember(constraint, ["inc", "dec"])
        diff = 1;
    else
        diff = 0;
    end

    v1 = zeros((nr_splines(1)-diff)*nr_splines(2), 1);
    for i=1:nr_splines(2)
        for j=1:nr_splines(1)-1
            v1(j+(i-1)*(nr_splines(1)-1)) = coef(j+(i-1)*nr_splines(1) + 1) - coef(j+(i-1)*nr_splines(1));
        end
    end

    if constraint == "inc"
        v1 = v1 < 0;
    elseif constraint == "dec"
        v1 = v1 > 0;
    elseif constraint == "none"
        v1 = zeros(size(v1));
    end

end

