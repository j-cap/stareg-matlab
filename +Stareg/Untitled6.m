%%  TEST PREDICTION VS FAST_PREDICTION
% data generation
n_data = 100000;
Xtrain = rand(n_data,2);
ytrain = 4*Xtrain(:,1).^2 - Xtrain(:,2) + randn(n_data, 1)*0.1;
% model description
d = {["s(1)", 100, "inc", 3000, "e"]; 
     ["t(1,2)", "12,20", "inc,dec", "2000,2000", "e,e"]};
%%
% build and train model
[coef, basis_matrix, model, reduced_model] = Stareg.fit(d, Xtrain, ytrain);

%% SAVE BOTH MODELS AND COMPARE SPACE
save("model.mat", "model");
save("reduced_model.mat", "reduced_model");

%% TEST PREDICTION 

n_pred = 5;
Xpred = rand(n_pred, 2);
ypred_true = 4*Xpred(:,1).^2 - Xpred(:,2);
figure(); scatter3(Xpred(:,1), Xpred(:,2), ypred_true); hold on; grid
%%
c = coef.c10;
%% TIME CHECK
tic; s_pred = Stareg.predict(Xpred, model, coef); toc;
tic; s_fast_prediction = Stareg.storage_efficient_prediction(Xpred, reduced_model, c); toc;

%%
figure(); 
scatter3(Xpred(:,1), Xpred(:,2), ypred_true); hold on; grid
scatter3(Xpred(:,1), Xpred(:,2), s_pred, 'x');
scatter3(Xpred(:,1), Xpred(:,2), s_fast_prediction, 'xg');

