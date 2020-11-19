%% Get some data
x = linspace(-6,6,100)';
y = 35*sin(x) + x.^3 + randn(size(x))*0.2;
y = y ./ 10;
scatter(x,y)

nparam = 25;
%% PLS fit using FD smoothness matrx
[b, X, S, k] = PLS(x, y, 1, nparam, 0);

%% PLS fit using analytic smoothness matrix
[b_a, X_a, S_a, k_a] = PLS(x, y, 1, nparam, 1);

%% Plot both fits

figure(); hold on;
title("PLS Fits");
scatter(x, y, 'DisplayName', 'Data');
plot(x, X * b, 'DisplayName', 'FD Smoothenss Matrix');
plot(x, X_a * b_a, 'DisplayName', 'Analytic Smoothenss Matrix');
legend();