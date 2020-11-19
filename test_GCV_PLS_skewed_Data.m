%% Get some data
rng(1);
x = randn(100, 1)' * 2;
y = 35*sin(x) + x.^3 + randn(size(x))*5;
y = y ./ 10;
scatter(x,y)

nparam = 25;
%% PLS fit using FD smoothness matrx
[b, X, best_lam k] = GCV_PLS(x, y, nparam, 100, 1, 0);

%% PLS fit using analytic smoothness matrix
[b_a, X_a, best_lam_a, k_a] = GCV_PLS(x, y, nparam, 100, 1, 1);

%% Plot both fits
figure(); hold on;
title("PLS Fits");
scatter(x, y, 'DisplayName', 'Data');
scatter(x, X * b, 'x', 'DisplayName', 'FD Smoothenss Matrix');
scatter(x, X_a * b_a, 's', 'DisplayName', 'Analytic Smoothenss Matrix');
legend();