rng(2);
%%
disp("Test the increasing constraint");
x = linspace(0, 2*pi, 100)';
y = 2.5*sin(x) + 0.25*x.^2 + 0.4*randn(size(x));

[b_c, X_c] = PIRLS(x,y,25,"inc", 1, 0);

%%
disp("Test the peak constraint");
x = linspace(0, 2, 50)';
y = exp(-(x - 0.4).^2 / 0.1) + 0.1*randn(size(x));

[b_c, X_c] = PIRLS(x,y,10,"peak", 1, 0);
%%
