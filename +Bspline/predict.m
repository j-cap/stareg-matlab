function s = predict(Xpred,knots,coef,ll)
%%
% Evaluate the B-spline given the knot sequence in knots with the
% coefficients in coef for data X.
%
% Parameters:
% X : array      - Data points of shape (n_samples, :) to evaluate the 
%                  B-spline on.
% knots : array  - Knot sequence, length = (length(coef) + ll + 1).
% coef : array   - B-spline coefficients.
% ll : double    - Order of the B-spline basis functions.
%
% Returns
% -------
% s : array      - B-spline values at X.

arguments
    Xpred (:,:) double
    knots (:,:) double
    coef (:,1) double
    ll (1,:) double = 3;
end

    [r,c] = size(Xpred);
    if c == 1
        disp('Prediction for 1-D data');
        B = zeros(r, length(coef));
        for i=ll+1:length(knots)-1 
            B(:,i-ll) = Bspline.basisfunction(Xpred, knots, i, ll);
        end
    elseif c == 2
        disp('Prediction for 2-D data');
        B1 = zeros(r, length(knots(:,1))-1-ll(1));
        B2 = zeros(r, length(knots(:,2))-1-ll(2));
        B = zeros(r, length(coef));
        
        for i=ll(1)+1:length(knots(:,1))-1
            B1(:,i-ll(1)) = Bspline.basisfunction(Xpred(:,1), knots(:,1), i, ll(1));
        end
        for i=ll(2)+1:length(knots(:,2))-1
            B2(:,i-ll(2)) = Bspline.basisfunction(Xpred(:,2), knots(:,2), i, ll(2));
        end
        for i=1:r
            B(i,:) = kron(B2(i,:), B1(i,:));
        end
    end
    
    s = B*coef;
end

