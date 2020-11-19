%% Increasing and Decreasing Data
x = linspace(0,6,100)';
y = 20*sin(x)*2 + x.^3 + randn(size(x))*5 ;
y = y / 10;

%% test increasing constraint
[beta_PLS_c, X] = PIRLS(x, y, 15, "inc", 1, 0);

%% test decreasing constraint
[beta_PLS_c, X] = PIRLS(x, -y, 15, "dec", 1, 0);

%% Peak and Valley Data
x = linspace(-5, 5, 100)';
y = 2*exp(-(x-1).^2 ./ 10) + 0.2*x + 0.1*randn(size(x));

%% test peak constraint
[beta_PLS_c, X] = PIRLS(x, y, 15, "peak", 1, 0);

%% test valley constraint
[beta_PLS_c, X] = PIRLS(x, -y, 15, "valley", 1,0);
