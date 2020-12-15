%% Try to implement the B-spline derivative from Fahrmeir, Regression p.430
ndata = 1000;
xdata = [linspace(-6,0,ndata), rand(1,ndata/4)*10];
xdata = sort(xdata)';
y = xdata.^3 - xdata.^2  + 4*cos(xdata) + 2*randn(size(xdata));
dy = 3*xdata.^2 - 2*xdata - 4*sin(xdata);
dyy = 6.*xdata - 2 - 4*cos(xdata);
%y = y ./ 100;
%dy = dy ./ 100;
%dyy = dyy ./ 100;
nsplines = 25;
l = 4;
[X, knots] = Bspline.basismatrix(xdata, nsplines, l, "q");
b = X \ y;
[dS, dX] = bspline_1st_derivative(xdata, knots, b, l);
[ddS, dX, dcoef, h2] = bspline_2nd_derivative_fahrmeir(xdata, knots, b, l);

%% Plot fit, first-order and second-order derivative
figure()
tiledlayout(3,1); nexttile
title("Data and Fit")
scatter(xdata, y); hold on;
plot(xdata, X * b);

nexttile; title("First Derivative")
plot(xdata, dy, ':', 'LineWidth', 2, 'DisplayName', 'analytic deriv.'); hold on;
plot(xdata, dS, '--', 'LineWidth', 2, 'DisplayName', 'spine deriv.');
legend()

nexttile; title("Second Derivative");
plot(xdata, dyy, 'DisplayName', 'analytic 2nd deriv.'); hold on;
plot(xdata, ddS, '--', 'LineWidth', 2, 'DisplayName', 'spline 2nd deriv. fahrmeir')
ylim([-40,40])
legend()

%%



