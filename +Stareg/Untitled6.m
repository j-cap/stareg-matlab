%%  TEST PREDICTION VS FAST_PREDICTION
% data generation
n_data = 100000;
Xtrain = rand(n_data,2);
ytrain = 4*Xtrain(:,1).^2 - Xtrain(:,2) + randn(n_data, 1)*0.1;
% model description
d = {["s(1)", 100, "dec", 3000, "e"]; 
     ["t(1,2)", "12,20", "inc,dec", "2000,2000", "e,e"]};
%%
% build and train model
[coef_, basis_matrix, model, reduced_model, coef_list] = Stareg.fit(d, Xtrain, ytrain);

%% SAVE BOTH MODELS AND COMPARE SPACE
save("model.mat", "model");
save("reduced_model.mat", "reduced_model");

%% TEST PREDICTION 

n_pred = 5;
Xpred = rand(n_pred, 2);
ypred_true = 4*Xpred(:,1).^2 - Xpred(:,2);
figure(); scatter3(Xpred(:,1), Xpred(:,2), ypred_true); hold on; grid
%%
c = coef_.c10;
%% TIME CHECK
tic; s_pred = Stareg.predict(Xpred, model, coef_); toc;
tic; s_pred_2 = Stareg.predict(Xpred, reduced_model, coef_); toc;

tic; s_fast_prediction = Stareg.storage_efficient_prediction(Xpred, reduced_model, coef_); toc;

%%
figure(); 
scatter3(Xpred(:,1), Xpred(:,2), ypred_true); hold on; grid
scatter3(Xpred(:,1), Xpred(:,2), s_pred, 'x');
scatter3(Xpred(:,1), Xpred(:,2), s_fast_prediction, 'xg');

