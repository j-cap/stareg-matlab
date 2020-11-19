function e = mse(y, ypred)
%% returns the mean squared error of y and ypred
e = sum((y - ypred).^2) / length(y);
end