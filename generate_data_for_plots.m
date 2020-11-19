%%

for i=4:17
    plot(xdata, bspline(xdata, knots3, i, l-1), 'DisplayName', ['Spline-', num2str(i-3)])
    hold on;
    legend()
end

%%

mean = 75;
sd = 5;
skew = .3;
kurt = 2.7;
col = 1
trial = 1;
rating = [];
while length(rating)<100
r = pearsrnd(mean,sd,skew,kurt,trial,col);
if r > 0 & r < 100
rating = [rating ; r];
end
end
hist(rating)

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
%%
%% Generate data for the spline-types plot
x = linspace(0,1,1000)';

knots_e = 0:0.1:1;
knots_q = [0, 0.1, 0.2, 0.38, 0.46, 0.52, 0.56, 0.83, 0.88, 0.9, 1];
knots = knots_e';
kidx = 1;
B0 = bspline(x, knots, kidx, 0);
B1 = bspline(x, knots, kidx+1, 1);
B2 = bspline(x, knots, kidx+2, 2);
B3 = bspline(x, knots, kidx+3, 3);

%%
k_save = [knots; zeros(numel(B0)-length(knots), 1)];
T = table(x, B0, B1, B2, B3, k_save, 'VariableNames', {'x', 'B0', 'B1', 'B2', 'B3', 'knots'});
writetable(T, "spline-types-quantile.txt", 'Delimiter', ' ');
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

%% Generate data for B-spline basis plots
ndata = 1000;
xe = linspace(0., 1, ndata)';
xq = exp(3*xe);
xq = sort((xq - min(xq)) / (max(xq) - min(xq)));

[Be, ke] = bspline_basis(xe, 8, 3, "e");
[Bq, kq] = bspline_basis(xq, 8, 3, "q");
%%
k_save_e = [ke'; zeros(ndata-length(ke), 1)];
k_save_q = [kq'; zeros(ndata-length(kq), 1)];

%dlmwrite("b-spline-basis-equidistant.txt", full([x, Be, k_save_e]), 'delimiter', ' ');
%dlmwrite("b-spline-basis-quantile.txt", full([x, Bq, k_save_q]), 'delimiter', ' ');

%%
figure();
plot(xq, Bq); hold on;
scatter(kq, zeros(size(kq)), 'o');



%% test the derivatives
x = [linspace(0,3,20), linspace(3.001, 6, 80)];
x = linspace(0,6,1000);
x = x';
y = sin(x.^2) + randn(size(x))*0.01;
yt = sin(x.^2);
sorder = 5;
[B, k] = bspline_basis(x, 35, sorder, "q");
b = B \ yt;

dB = bspline_1st_derivative(x, k, b, sorder);
[dB2eilers, dXeilers] = bspline_2nd_derivative_eilers(x, k, b, sorder);
[dB2fahrm, dXfahrm] = bspline_2nd_derivative_fahrmeir(x, k, b, sorder);
dy = 2.*cos(x.^2).*x;
dy2 = -4.*sin(x.^2).*x.^2 + 2*cos(x.^2);

figure();
scatter(x, y); hold on;
plot(x, B * b, 'LineWidth', 2);
title("Data and Fit");
legend("Data", "B-spline");

figure(); 
plot(x, dB); hold on;
plot(x, dy);
title("First derivative");
legend("spline", "analyitic")

figure(); 
plot(x, dB2eilers, 'LineWidth', 2); hold on;
plot(x, dy2, 'LineWidth', 2);
plot(x, dB2fahrm, 'LineWidth', 2);
title("Second derivative");
legend("Eilers", "analytic","Fahrmeir");

%%
P = dXfahrm' * dXfahrm;
image(P, 'CDataMapping', 'scaled');

%% Test the TPS P-spline formulation in Fahrmeir, p.508
d1 = 5; d2 = 5;
I1 = eye(d1); I2 = eye(d2);
D1 = full(smoothness_matrix(d1, 2, 0));
D2 = full(smoothness_matrix(d2, 2, 0));

K1 = kron(I2, D1)' * kron(I2, D1);
K2 = kron(D2, I1)' * kron(D2, I1);

b = [1:1:25]';
b = b.^2;
br = reshape(b, 5,5);
% calculate penalty term in matrix notation
b' * (K1 + K2) * b
%%
bb = 0;
% calculate penalty term in sum notation
for i=1:d1
    for j=3:d2
        bb = bb + (br(i,j) - 2*br(i, j-1) + br(i, j-2))^2;
    end
end

for i=3:d1
    for j=1:d2
        bb = bb + (br(i,j) - 2*br(i-1, j) + br(i-2, j))^2;
    end
end
bb

%%
d1 = 10; d2 = 5;
I1 = eye(d1); I2 = eye(d2);
D12 = smoothness_matrix(d1, 2, 0); D22 = smoothness_matrix(d2, 2, 0);

K1 = kron(I2, D12)' * kron(I2, D12);
K2 = kron(D22, I1)' * kron(D22, I1);

b' * (K1 + K2) * b

%%
bb = 0;
br = reshape(b, d1, d2);
% calculate penalty term in sum notation
for i=1:d1
    for j=3:d2
        bb = bb + (br(i,j) - 2*br(i, j-1) + br(i, j-2))^2;
    end
end

for i=3:d1
    for j=1:d2
        bb = bb + (br(i,j) - 2*br(i-1, j) + br(i-2, j))^2;
    end
end
bb

%%
for i=1:11
surf(xg, yg,reshape(X(:,i)*b(i), nd, nd)); hold on; pause(1);
end



