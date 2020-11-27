%% Try to implement the B-spline derivative from Fahrmeir, Regression p.430
ndata = 100;
xdata = linspace(-6,6,ndata);
%xdata = randn(1, ndata)*10 + randn(1, ndata)*0.2;
xdata = sort(xdata)';
y = xdata.^3 - xdata.^2  + 2*randn(size(xdata)) + 4*cos(xdata);
dy = 3*xdata.^2 - 2*xdata + -4*sin(xdata);
dyy = 6.*xdata - 2;
nsplines = 55;
l = 3;
[X, knots] = bspline_basis(xdata, nsplines, l, "q");
b = X \ y;
%%
scatter(xdata, y)
%% Plot data, fit and individual splines
figure(); hold on;
scatter(xdata(1:2:end), y(1:2:end));
plot(xdata, X * b);
legend("data", "fit");
for i=1:nsplines
    plot(xdata, X(:,i)*b(i), 'DisplayName', ['Spline-', num2str(i)]);
end

%% Plot analytic and 1st spline derivative for cubic splines
[X3, knots3] = bspline_basis(xdata, nsplines, 3, "q");
b3 = X3 \ y;
dyfit3_FAHRM = bspline_1st_derivative(xdata, knots3, b3, 3);
dyfit3_HOFNER = bspline_1st_derivative_HOFNER(xdata, knots3, b3, 3);
figure(); title("Fit and first derivative")
plot(xdata, dy); hold on;
plot(xdata, dyfit3_FAHRM, 'x');
plot(xdata, dyfit3_HOFNER);
scatter(xdata(1:2:end), y(1:2:end));
plot(xdata, X3 * b3);
legend("analytisch deriv", "fahrmeir deriv", "hofner deriv", "data", "fit");


%% Plot analytic and 1nd spline derivative for quadratic splines
[X2, knots2] = bspline_basis(xdata, nsplines, 2, "e");
b2 = X2 \ y;
dyfit2 = bspline_1st_derivative(xdata, knots2, b2, 2);
figure(); title("Fit and first derivative")
plot(xdata, dy); hold on;
plot(xdata, dyfit2);
scatter(xdata(1:20:end), y(1:20:end));
plot(xdata, X * b);
legend("analytisch deriv", "spline deriv", "data", "fit");

%% 
ndata = 100;
xdata = linspace(0,2,ndata);
y = sin(xdata)';
nsplines = 21;
l = 3;
[X3, k3] = bspline_basis(xdata, nsplines, l, "e");
b = X3 \ y;
length(k3) - nsplines - 1
%%
DX = bspline_1st_derivative(xdata, k3, b, l);

%% calculate the first derivative as in smoothness_matrix (analytic part)
xx = xdata;
dS = zeros(length(xx), nsplines);
for i=l+1:length(k3)-1
    left = bspline(xx, k3, i-1, l-1);
    right = bspline(xx, k3, i, l-1);
    dS(:,i-l) = l * (left - right) / (mean(diff(k3))*l);
end 
plot(dS)
%% calculate the second derivative using the same scheme as above
xx = xdata;
ddS = zeros(length(xx), nsplines);
for i=l+1:length(k3)-1
    left = bspline(xx, k3, i-1, l-2);
    right = bspline(xx, k3, i, l-2);
    ddS(:,i-l) = l * (left - right) / (mean(diff(k3))*l);
end 
plot(ddS)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test GCV for PLS and PLS_analytic
x = linspace(0,6,1000)';
y = 40*sin(x)*2 + x.^3 + randn(size(x))*5 ;
y = y / 10;
scatter(x,y)
%% Bspline basis and smoothness matrices for FD and analytic derivative
[X, k] = bspline_basis(x, 25, 2, "e");
S_fd = smoothness_matrix(25, 2, 0);
S_ana = smoothness_matrix(25, 2, k);
%% GCV for both types
[b_gcv, X_gcv, best_lam] = GCV_PLS(x, y, 25, 100, 1, 0);
[b_gcv_ana, X_gcv_ana, best_lam_ana] = GCV_PLS(x, y, 25, 100, 1, 1);
%% Plotting
figure();
ip = 20;
scatter(x(1:ip:end), y(1:ip:end)); hold on;
plot(x, X_gcv * b_gcv, 'LineWidth', 2);
plot(x, X_gcv_ana * b_gcv_ana, 'LineWidth', 2);
legend('Data', 'FD Smoothing', 'Analytic Smoothing');
disp(["Optimal Lambda FD = ", best_lam]);
disp(["Optimal Lambda Analytic= ", best_lam_ana]);

%%
disp(["MSE FD Smoothing = ", mse(y, X_gcv * b_gcv)]);
disp(["MSE Analytic Smoothing = ", mse(y, X_gcv_ana * b_gcv_ana)]);



