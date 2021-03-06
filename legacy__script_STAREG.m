%% get some data
x = linspace(-3,3,100);
[xg, yg] = meshgrid(x);

data = [xg(:) yg(:)];

f = @(x,y) 2*exp(-0.1*(x - 0.5).^2 - 2*(y-0.5).^2) + randn(size(x))*0.1;
z = f(data(:,1), data(:,2));

scatter3(data(:,1), data(:,2), z);
%% get TPS basis
nparam1 = 25;
nparam2 = 25;
[X, k1, k2] = tpspline_basis(data, nparam1, nparam2, 2);
%% get Smoothness Matrix for TPS
S = smoothness_matrix_TPS(nparam1, nparam2, 2);
%% solve PLS for TPS
beta_LS = (X'*X) \ X' * z;
beta_PLS = ((X'*X) + 0.1 * (S'*S)) \ X' * z;
%% test PLS_tps.m
[beta_PLS_func, XX, SS] = PLS_tps(data, z, 0.1, nparam1, nparam2);
%% Plot
figure(); 
scatter3(data(:,1), data(:,2), z); hold on; 
scatter3(data(:,1), data(:,2), X * beta_LS);
scatter3(data(:,1), data(:,2), X * beta_PLS);
scatter3(data(:,1), data(:,2), X * beta_PLS_func);
legend('Data','LS','PLS', 'PLS_func');

%%
[beta_PLS_GCV, X, best_lam] = GCV_PLS(data, z, [12, 20], 100, 0, 0);

%% STRUCTURED ADDITIVE REGRESSION
n_data = 1000;
X = linspace(0,1,150)';
XX = [rand(n_data,1), rand(n_data,1)];
% XX = [linspace(0,1,n_data)+5; linspace(0,1,n_data)]';
ytrue = 4*X.^2 + sin(X.*6);
y = ytrue + randn(size(X))*0.1;
yy_true = 0.5*XX(:,1).^2 - 2*cos(XX(:,2)*6) + XX(:,1).*XX(:,2); 
yy = yy_true + randn(size(XX(:,1)))*0.4;
%% 1D
figure();
plot(X, ytrue, "DisplayName", "True Function"); hold on;
scatter(X, y, "DisplayName", "Data");
scatter(X, BB*cc, "DisplayName", "Fit");
legend; grid;
%% 2D
figure()
scatter3(XX(:,1), XX(:,2), yy, 'DisplayName','Data'); hold on;
scatter3(XX(:,1), XX(:,2), yy_true, 'DisplayName','True Function'); hold on;
scatter3(XX(:,1), XX(:,2), BB*cc, 'DisplayName', 'Fit');
legend()
%% create the cell array description of the model
d = {["s(1)", "inc", 80, 6000, "e"];
     ["s(2)", "peak", 80, 6000, "e"];
     ["t(1,2)", "inc,none", "12,8", "6000,6000", "e,e"]};
%%      
[cc, BB] = Stareg.fit(d, XX, yy);
%% 1D
figure();
plot(X, ytrue, "DisplayName", "True Function"); hold on;
scatter(X, y, "DisplayName", "Data");
scatter(X, BB*cc, "DisplayName", "Fit");
legend; grid;
%%

[w, wc] = Utils.check_constraint_full(coef_pls, m, B, yy);
%%
iterIdx = 0;
BtB = B'*B;
Bty = B'*yy;

W_compare = W;
W_old = zeros(size(W));


while ~all(W_old == W_compare)
    W_old = W_compare;
    coef_pls = (BtB + S + K) \ Bty;
    
    [w, W_compare] = Utils.check_constraint_full(coef_pls, m, B, yy);
    K = Utils.create_constraint_matrix(m, w);
    
    iterIdx = iterIdx + 1;
    iterIdx
    if iterIdx > 5
        disp("Stop the count")
        break
    end
end
%%
image(KKK, 'CDataMapping','scaled'); colorbar




%%
figure();
plot(X, ytrue, "DisplayName", "True Function"); hold on;
scatter(X, y, "DisplayName", "Data");
scatter(X, B*coef, "DisplayName", "Pspline");
legend; grid;
%%
for i=1:numel(d)
    c = d{i};
    model(i).type = c{1};
    model(i).constraint = c{2}; 
    if model(i).type(1) == 's'
        nsplines = [str2num(c(3))];
        j = 4;
    elseif model(i).type(1) == 't'
        nsplines = [str2num(c(3)), str2num(c(4))];
        j = 5;
    end
    model(i).nsplines = nsplines;
    model(i).lam_s = str2double(c(j));
    model(i).lam_c = str2double(c(j+1));
    model(i).knot_type = c{j+2};
end

%%


