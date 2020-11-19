%% test for bspline_basis.m

xdata = linspace(0, 2*pi, 200);
nsplines = 10;
[X0, k0] = bspline_basis(xdata, nsplines, 0);
[X1, k1] = bspline_basis(xdata, nsplines, 1);
[X2, k2] = bspline_basis(xdata, nsplines, 2);
[X3, k3] = bspline_basis(xdata, nsplines, 3);

figure();
%% cubic splines
for i=1:nnsplines
    plot(xdata, X3(:,i)); hold on;
end
%% quadratic splines
for i=1:nnsplines
    plot(xdata, X2(:,i)); hold on;
end
%% linear splines 
for i=1:nnsplines
    plot(xdata, X1(:,i)); hold on;
end
%% constant splines
for i=1:nnsplines
    plot(xdata, X0(:,i)); hold on;
end

