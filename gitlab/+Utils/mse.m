function e = mse(y, ypred)
%% returns the mean squared error of y and ypred
arguments
    y (:,1) double
    ypred (:,1) double
end
e = sum((y - ypred).^2) / length(y);
end