


%% generate some data
ndata = 2000;
x = rand(ndata, 1);
y = rand(ndata, 1);
f = @(x,y) 2*exp(-(x - 0.2).^2 / 0.01) + 2.*y + 0.1*randn(size(x));
z = f(x,y);

ngrid=100;
xl = linspace(0,1,ngrid);
[xg,yg] = meshgrid(xl,xl); zg = f(xg, yg);

figure();
scatter3(x, y, z, 400, 'r'); hold on;
surf(xg,yg,zg);
%% tensorproductsplines for standard data
nsplines = 12;
Bx = bspline_basis(x, nsplines, 3, "e");
By = bspline_basis(y, nsplines+2, 3, "e");

X = zeros(ndata, nsplines*(nsplines+2));
for i=1:ndata
    X(i,:) = kron(By(i,:), Bx(i,:));
end
X = sparse(X);
b = X \ z;
figure(); 
scatter3(x,y,z); hold on; 
scatter3(x,y,X * b);
legend('Data', 'TPS Fit');


%% tensorproductsplines for gridded data
%xtg = unique(xg);
%ytg = unique(yg);

nsplines = 12;
Bxg = bspline_basis(xg(:), nsplines, 2);
Byg = bspline_basis(yg(:), nsplines, 2);

Xg = zeros(ngrid*ngrid, nsplines*nsplines);
for i=1:ngrid*ngrid
    Xg(i,:) = kron(Bxg(i,:), Byg(i,:));
end
Xg = sparse(Xg);

bg = Xg \zg(:);
figure(); 
scatter3(xg(:),yg(:),zg(:)); hold on; 
scatter3(xg(:),yg(:),Xg * bg);
legend('Data', 'TPS Fit');

%% Plot Tensorproductsplinebasis
figure();
for i=1:200:205
    surf(xg, yg, reshape(Xg(:, i), size(xg))); hold on;
end

%%
nd = 50;
x = linspace(0,1,nd)';
y = linspace(0,10,nd)';
[xg, yg] = meshgrid(x,y);
zg = sin(yg) + 3*xg.^2;
[X, k1, k2] = tpspline_basis([xg(:), yg(:)], 10, 5, 3);

b = X \ zg(:);

%%
figure();
scatter3(xg(:),yg(:), zg(:)); hold on;
scatter3(xg(:),yg(:), X*b, "rx");
idx = 1;
surf(xg, yg, reshape(X(:,idx), nd, nd)*b(idx)); hold on
legend('data', 'fit', 'individual spline');
%%
for i=1:10
    surf(xg, yg, reshape(X(:,i)*b(i), nd, nd)); hold on
end


