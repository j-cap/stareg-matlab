%%
xr = randn(1000,1);
y = sin(xr);


%% equidistant knots
l = 3;
[Xe, ke] = bspline_basis(xr, 10, l, "e");
be = Xe \ y;

figure(); hold on;
for i=1:numel(be)
    scatter(xr, Xe(:,i), 'x');
end
for i=1:numel(ke)
    lw = 1;
    if i == l+1 | i == length(ke) - l lw = 2; end
    xline(ke(i), 'LineWidth', lw);
    scatter(ke(i), 0, 'ok');
end
figure();
scatter(xr, y, 'x', 'DisplayName', 'Data'); hold on;
scatter(xr, Xe * be, 'DisplayName', 'Fit');
legend();

%% quantile knots
l = 3;
[Xq, kq] = bspline_basis(xr, 50, 3, "q");
bq = Xq \ y;

figure(); hold on;
for i=1:numel(bq)
    scatter(xr, Xq(:,i), 'x');
end
for i=1:numel(kq)
    lw = 1;
    if i == l+1 | i == length(kq) - l lw = 2; end
    xline(kq(i), 'LineWidth', lw);
    scatter(kq(i), 0, 'ok');
end

%%
figure();
scatter(xr, y, 'x', 'DisplayName', 'Data'); hold on;
scatter(xr, Xq * bq, 'DisplayName', 'Fit');
legend();


