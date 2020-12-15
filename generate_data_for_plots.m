%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate data for the Bias-Variance Plot
x = linspace(0.1,2,250)';
y_bias = 5 ./ (2*x) + 2;
y_var = 4*x.^3 + 2;
y_both = y_bias + y_var;
[optim, o_idx] = min(y_both);
%%
%T = table(x, y_bias, y_var, y_bias+y_var, ones(size(x))*x(o_idx), 'VariableNames', {'x', 'bias', 'variance', 'both', 'optimum'});
%writetable(T, "bias-variance-decomposition.txt", 'Delimiter', ' ');
plot(x, y_bias); hold on;
plot(x, y_var);
plot(x, y_bias + y_var)
legend('bias', 'variance', 'both');
xline(x(o_idx));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate data for the spline-types plot
x = linspace(0,1,1000)';

knots_e = 0:0.1:1;
knots_q = [0, 0.1, 0.2, 0.38, 0.46, 0.52, 0.56, 0.83, 0.88, 0.9, 1];
knots = knots_e';
kidx = 1;
B0 = Bspline.basisfunction(x, knots, kidx, 0);
B1 = Bspline.basisfunction(x, knots, kidx+1, 1);
B2 = Bspline.basisfunction(x, knots, kidx+2, 2);
B3 = Bspline.basisfunction(x, knots, kidx+3, 3);

%%
%k_save = [knots; zeros(numel(B0)-length(knots), 1)];
%T = table(x, B0, B1, B2, B3, k_save, 'VariableNames', {'x', 'B0', 'B1', 'B2', 'B3', 'knots'});
%writetable(T, "spline-types-quantile.txt", 'Delimiter', ' ');
%%
figure();
plot(x, B0); hold on;
plot(x, B1); 
plot(x, B2);
plot(x, B3);
legend('B0', 'B1', 'B2', 'B3');
scatter(knots, zeros(size(knots)));
for k=1:numel(knots)
    xline(knots(k));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate data for B-spline basis plots
ndata = 1000;
xe = linspace(0., 1, ndata)';
xq = exp(3*xe);
xq = sort((xq - min(xq)) / (max(xq) - min(xq)));

[Be, ke] = Bspline.basismatrix(xe, 8, 3, "e");
[Bq, kq] = Bspline.basismatrix(xq, 8, 3, "q");
%%
k_save_e = [ke'; zeros(ndata-length(ke), 1)];
k_save_q = [kq'; zeros(ndata-length(kq), 1)];

%dlmwrite("b-spline-basis-equidistant.txt", full([x, Be, k_save_e]), 'delimiter', ' ');
%dlmwrite("b-spline-basis-quantile.txt", full([x, Bq, k_save_q]), 'delimiter', ' ');

%%
figure();
plot(xe, Be); hold on;
scatter(kq, zeros(size(kq)), 'o');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
