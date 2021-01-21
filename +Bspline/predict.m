function s = predict(Xpred,knots,coef,order)
%%
% Evaluate the B-spline given the knot sequence in knots with the
% coefficients in coef for data X.
%
% Parameters:
% X : array      - Data points of shape (n_samples, :) to evaluate the 
%                  B-spline on.
% knots : array  - Knot sequence, length = (length(coef) + order + 1).
% coef : array   - B-spline coefficients.
% order : double - Order of the B-spline basis functions.
%
% Returns
% -------
% s : array      - B-spline values at X.

arguments
    Xpred (:,:) double
    knots (:,:) struct
    coef (:,1) double
    order (1,:) double = 3;
end

    [r,c] = size(Xpred);
    if c == 1
        disp('Prediction for 1-D data');
        B = zeros(r, length(coef));
        for i=order+1:length(knots.k1)-1 
            B(:,i-order) = Bspline.basisfunction(Xpred, knots.k1, i, order);
        end
    elseif c == 2
        disp('Prediction for 2-D data');
        B1 = zeros(r, length(knots.k1)-1-order(1));
        B2 = zeros(r, length(knots.k2)-1-order(2));
        B = zeros(r, length(coef));
        
        for i=order(1)+1:length(knots.k1)-1
            B1(:,i-order(1)) = Bspline.basisfunction(Xpred(:,1), knots.k1, i, order(1));
        end
        for i=order(2)+1:length(knots.k2)-1
            B2(:,i-order(2)) = Bspline.basisfunction(Xpred(:,2), knots.k2, i, order(2));
        end
        for i=1:r
            B(i,:) = kron(B2(i,:), B1(i,:));
        end
    end
    
    s = B*coef;
end

