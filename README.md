# stareg-matlab: STructured  Additive REGression in Matlab


## Usage

```Matlab

rng(2); % set random seed
n_data = 500; 
X = [rand(n_data,1), rand(n_data,2)]; % create input data
y = exp(-(X(:,1)-0.4).^2 ./ 0.01) + X(:,2).^2 + X(:,1).*X(:,2) + randn(n_data,1)*0.2; % create target data


d = {["s(1)", "inc", 100, 3000, "e"]; 
     ["t(1,2)", "none,inc", "12,20", "2000,2000", "e,q"]}; % create model description
     
[coef, basis_matrix, model] = Stareg.fit(d, X, y); % fit the model
```

## Features
- Structured Additive Regression
- B-Splines, P-Splines and shape-constraint P-splines
- Incorporate prior knowledge using constraints:
  - Monotonicity constraints: increasing, decreasing
  - Shape constraints: convex, concave
  - Peak and Valley constraints
- Interpretable Machine Learning
- Uni- and multivariate regression


