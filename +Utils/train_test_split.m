function [Xtrain, ytrain, Xtest, ytest] = train_test_split(X,y,test_size)
%%
% Perform a train-test split on the given data X, y.
%
% Parameters:
% -----------
% X : array             - Input data of shape (n_samples, n_dim).
% y : array             - Output data of shape (n_samples, 1).
% test_size : double    - Size of the test set.
%
% Returns:
% --------
% Xtrain : array 
% ytrain : array 
% Xtest : array
% ytest : array
arguments
    X (:,:) double
    y (:,1) double
    test_size (1,1) double = 0.3
end
rng(2);
N = size(X,1);
idx = randperm(N);
Xtest = X(idx(1:round(N*test_size)),:);
ytest = y(idx(1:round(N*test_size)),:);

Xtrain = X(idx(round(N*test_size)+1:end),:);
ytrain = y(idx(round(N*test_size)+1:end),:);

end

