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

d = {["s(1)", "inc", 20, 0.1, 300, "equidistant"]; 
     ["t(1,2)", "none", [12,23], 0.2, 2000, "quantile"];
     ["s(2)", "inc", 20, 0.1, 300, "equidistant"]; 
     ["s(3)", "dec", 20, 0.1, 300, "equidistant"]; 
     ["s(4)", "inc", 0, -0.1, 300, "equidistant"]; };
     
model = {};
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

