%%
rng(2); % set random seed
n_data = 500; 
X = [rand(n_data,1), rand(n_data,2)]; % create input data
y = exp(-(X(:,1)-0.4).^2 ./ 0.01) + X(:,2).^2 + X(:,1).*X(:,2) + randn(n_data,1)*0.2; % create target data

Xtest = [rand(100,1), rand(100,1)]; % generate test data
ytest = exp(-(Xtest(:,1)-0.4).^2 ./ 0.01) + Xtest(:,2).^2 + Xtest(:,1).*Xtest(:,2)

d = {["s(1)", "peak", 100, 3000, "e"]; 
     ["t(1,2)", "none,inc", "12,20", "2000,2000", "e,q"]}; % create model description
     
[coef, basis_matrix, model] = Stareg.fit(d, X, y); % fit the model
ypred = Stareg.predict(Xtest, model, coef);

fprintf("MSE on test data = %4.4f \n", Utils.mse(ytest, ypred));

%%
rng(2);

%% 2d example
ndata = 200;
X = [rand(ndata,1), rand(ndata,1)];
f = @(x1,x2) x1.^2 + sqrt(x2) - x1.*x2*0.5;
y = f(X(:,1), X(:,2)) + randn(ndata,1)*0.05;
%%
npred = 25;
Xpred = [rand(npred,1), rand(npred,1)];
ypred = f(Xpred(:,1), Xpred(:,2));
scatter3(X(:,1), X(:,2), y);
%%
d = {%["s(1)", "inc", 100, 3000, "e"], 
    %["s(2)", "none", 100, 3000, "e"]; 
    ["t(1,2)", "inc,none", "12,20", "2000,2000", "e,q"]};
[coef, BB, model] = Stareg.fit(d, X, y);
%%

y__PREDICTION = Stareg.predict(Xpred, model, coef);
Utils.mse(y__PREDICTION, ypred)

%%
%Xpred = xpred;
i = 1;
fn = fieldnames(model);
nprediction = length(Xpred);
B = [];

for i=1:numel(fn)

    type_ = model.(fn{i}).type
    nr_coef = model.(fn{i}).nr_splines; % get coefficients of the submodel
    knots = model.(fn{i}).knots;     % get knots of the submodel
    order = model.(fn{i}).spline_order; % get spline order of the submodel
   
    if type_.startsWith("s") % 
        dim = str2double(type_{1}(3));
        %sub_coef = coef__(i:i+nr_coef-1); 
        %sub_prediction = Bspline.predict(Xpred(:,dim), knots, sub_coef, order);
        disp('Prediction for 1-D data');
        basis = zeros(nprediction, nr_coef);
        for s_idx=order.so1+1:length(knots.k1)-1 
            basis(:,s_idx-order.so1) = Bspline.basisfunction(Xpred(:,dim), knots.k1, s_idx, order.so1);
        end
        
    elseif type_.startsWith("t") % check the constraint for tensor-product B-splines
        dim1 = str2double(type_{1}(3));
        dim2 = str2double(type_{1}(5));
        %sub_coef = coef__(i:i+prod(nr_coef)-1);
        %sub_prediction = Bspline.predict(Xpred(:,[dim1,dim2]), knots, sub_coef, order);
        
        disp('Prediction for 2-D data');
        B1 = zeros(nprediction, length(knots.k1)-1-order.so1); %basis matrix for direction 1
        B2 = zeros(nprediction, length(knots.k2)-1-order.so2); %basis matrix for direction 2
        basis = zeros(nprediction, prod(nr_coef));
        
        for s_idx1=order.so1+1:length(knots.k1)-1
            B1(:,s_idx1-order.so1) = Bspline.basisfunction(Xpred(:,dim1), knots.k1, s_idx1, order.so1);
        end
        for s_idx2=order.so2+1:length(knots.k2)-1
            B2(:,s_idx2-order.so2) = Bspline.basisfunction(Xpred(:,dim2), knots.k2, s_idx2, order.so2);
        end
        for data_idx=1:nprediction
            basis(data_idx,:) = kron(B2(data_idx,:), B1(data_idx,:));
        end
    else
        disp("Only B-splines (s) and tensor-product B-splines (t) available");
    end
    B = [B, basis];
    i = i + prod(nr_coef); 
end

s = B*coef__;
%%
%%
figure()
scatter3(X(:,1), X(:,2), y, 'DisplayName', 'Data'); hold on
scatter3(X(:,1), X(:,2), BB*coef.("c"+string(length(fieldnames(coef))-1)), 'DisplayName', 'Data'); 
%%
figure()
scatter3(Xpred(:,1), Xpred(:,2), ypred, 'DisplayName', 'Data'); hold on; 
scatter3(Xpred(:,1), Xpred(:,2), s, 'DisplayName', 'Fit'); legend()

%%
figure()
scatter(x,y); hold on; 
scatter(x, basis*coef.c0, 'DisplayName', 'PLS');
scatter(x, basis*coef.c4, 'DisplayName', 'Fit');
%%
figure()
scatter(xpred, ypred, 'x'); hold on; scatter(xpred, s, 'DisplayName', 'Fit');



%%
%% 1d example
x = rand(ndata,1);
y = 2*exp(-(x - 0.4).^2 / 0.02) + randn(ndata,1)*0.1;

xpred = linspace(0,1,100)';
ypred = 2*exp(-(xpred - 0.4).^2 / 0.02);

scatter(x,y); hold on; scatter(xpred, ypred, 'x');
d = {["s(1)", "peak", 100, 3000, "e"]}; 
     %["t(1,2)", "dec,none", "12,20", "2000,2000", "e,q"]};
%% 
[coef, basis, model] = Stareg.fit(d, x, y);
%%
