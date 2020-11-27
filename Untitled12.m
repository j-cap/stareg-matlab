
%%

x = linspace(-3,3,25);
y = linspace(0,6,25);
[xg, yg] = meshgrid(x,y);
zg = 5*exp(-xg.^2 ./ 10 - (yg-3.4).^2) + randn(size(xg))*0.1;

figure()
surf(xg,yg,zg, 'FaceAlpha', 0.1); %alpha = 0.1;
hold on;

X = [xg(:), yg(:)];
nr_splines = [7,10];
best_lam = Bspline.calc_GCV(X, zg(:), 10, 0, nr_splines);
[coef, T, k1, k2] = Bspline.fit_tp_Pspline(X, zg(:), best_lam, nr_splines);

scatter3(X(:,1), X(:,2), T*coef, 'rx');

%%
j = 14;
for i=j*nr_splines(1):(j+1)*nr_splines(1)
    surf(xg,yg,coef(i)*reshape(T(:,i), 25, 25)); hold on;
end