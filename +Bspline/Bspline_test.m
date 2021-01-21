function tests = Bspline_test
%%
% Main Function 
tests = functiontests(localfunctions);
end

function test_basisfunction(testCase)
    n_data = 25;
    x = linspace(0,1,n_data)';
    knots = linspace(0,1,5);
    b4_3 = Bspline.basisfunction(x, knots, 4, 3);
    assert(all(size(x) == size(b4_3)));
end

function test_basismatrix_equidistant(testCase)
    n_data = 100;
    x = linspace(0,1,n_data)';
    nr_splines = 10;
    order = 3;
    [basis_matrix, knots] = Bspline.basismatrix(x, nr_splines, order, "e");
    assert(all(size(basis_matrix) == [n_data, nr_splines]));
    assert(length(knots.k1) == nr_splines+order+1);
    assert(Utils.iscloseall(sum(full(basis_matrix), 2), ones(n_data,1)));
    Utils.iscloseall(mean(diff(knots.k1)), knots.k1(2)-knots.k1(1))
    Utils.iscloseall(mean(diff(knots.k1)), knots.k1(end-2)-knots.k1(end-1))
end

function test_basismatrix_quantile(testCase)
    n_data = 100;
    x = linspace(0,1,n_data)';
    nr_splines = 10;
    order = 3;
    [basis_matrix, knots] = Bspline.basismatrix(x, nr_splines, order, "q");
    assert(all(size(basis_matrix) == [n_data, nr_splines]));
    assert(length(knots.k1) == nr_splines+order+1);
    assert(Utils.iscloseall(sum(basis_matrix, 2), ones(n_data, 1)));
end

function test_tensorproduct_basismatrix(testCase)

    n_data = 100;
    x1 = rand(n_data,1);
    x2 = linspace(0,1,n_data)';
    X = [x1, x2];
    orders = [3,3];
    knot_types = ["e", "q"];
    nr_splines = [12, 8];
    [basis_matrix, knots] = Bspline.tensorproduct_basismatrix(X, nr_splines, orders, knot_types);
    
    assert(length(knots.k1) == nr_splines(1)+orders(1)+1);
    assert(length(knots.k2) == nr_splines(2)+orders(2)+1);
    assert(all(size(basis_matrix) == [n_data, prod(nr_splines)]));
    assert(Utils.iscloseall(sum(basis_matrix, 2), ones(n_data, 1)));
    
end

function test_fit_Bspline(testCase)
    rng(2); % set seed for random number generator
    n_data = 100;
    x = linspace(0,1,n_data)';
    y = sin(x*6) + randn(n_data,1)*0.2;
    nr_splines = 50;
    order = 3;
    knot_type = "e";
    [coef, basis_matrix, knots] = Bspline.fit(x,y,nr_splines, order, knot_type);
    
    assert(length(coef) == nr_splines);
    assert(all(size(basis_matrix) == [n_data, nr_splines]));
    assert(Utils.iscloseall(sum(basis_matrix, 2), ones(n_data, 1)));
    assert(length(knots.k1) == nr_splines+order+1);

    figure();
    title("Bspline Test"); hold on;
    scatter(x,y, 'DisplayName', 'Data');
    plot(x, basis_matrix*coef, 'DisplayName', 'Fit'); 
    legend(); grid();
    
end

function test_fit_tensorproduct_Bspline(testCase)
    rng(2); % set seed for random number generator
    n_data = 100;
    x1 = rand(100,1);
    x2 = rand(100,1);
    [x1g, x2g] = meshgrid(x1,x2);
    X = [x1g(:), x2g(:)];
    y = sin(X(:,1)*6) + 2*X(:,2).^2 + randn(n_data^2,1)*0.2;
    
    nr_splines = [12,8];
    orders = [3,3];
    knot_types = ["e", "e"];
    
    [coef, basis_matrix, knots] = Bspline.fit(X,y, nr_splines, orders, knot_types);
    
    assert(length(coef) == prod(nr_splines));
    assert(all(size(basis_matrix) == [n_data^2, prod(nr_splines)]));
    assert(Utils.iscloseall(sum(basis_matrix, 2), ones(n_data^2, 1)));
    assert(length(knots.k1) == nr_splines(1)+orders(1)+1);
    assert(length(knots.k2) == nr_splines(2)+orders(2)+1);

    figure(); 
    scatter3(X(:,1), X(:,2), y); hold on; 
    scatter3(X(:,1), X(:,2), basis_matrix*coef);
    legend("Data", "Fit");  
    title("Tensorproduct Bspline Test");
    
end

function test_fit_Pspline(testCase)
    rng(2); % set seed for random number generator
    n_data = 100;
    x = linspace(0,1,n_data)';
    y = sin(x*6) + randn(n_data,1)*0.2;
    nr_splines = 50;
    order = 3;
    knot_type = "e";
    lam = 1;
    [coef, basis_matrix, knots] = Bspline.fit_Pspline(x,y,lam,nr_splines, order, knot_type);
    
    assert(length(coef) == nr_splines);
    assert(all(size(basis_matrix) == [n_data, nr_splines]));
    assert(Utils.iscloseall(sum(basis_matrix, 2), ones(n_data, 1)));
    
    figure();
    title("Pspline Test"); hold on;
    scatter(x,y, 'DisplayName', 'Data');
    plot(x, basis_matrix*coef, 'DisplayName', 'Fit'); 
    legend(); grid();
    
    % test if Bspline.fit_Pspline() with lam == 0 produces the same as
    % Bspline.fit()
    [coef_bs, basis_matrix_bs, knots_bs] = Bspline.fit(x,y, nr_splines, order, knot_type);
    [coef_ps0, basis_matrix_ps0, knots_ps0] = Bspline.fit_Pspline(x,y,0,nr_splines, order, knot_type);
    assert(Utils.iscloseall(coef_bs, coef_ps0));
    
end


function test_fit_tensorproduct_Pspline(testCase)
    rng(2); % set seed for random number generator
    n_data = 100;
    x1 = rand(100,1);
    x2 = rand(100,1);
    [x1g, x2g] = meshgrid(x1,x2);
    X = [x1g(:), x2g(:)];
    y = sin(X(:,1)*6) + 2*X(:,2).^2 + randn(n_data^2,1)*0.2;
    
    nr_splines = [12,8];
    orders = [3,3];
    knot_types = ["e", "e"];
    lam = 1;
    
    [coef, basis_matrix, knots] = Bspline.fit_Pspline(X,y, lam, nr_splines, orders, knot_types);
    
    assert(length(coef) == prod(nr_splines));
    assert(all(size(basis_matrix) == [n_data^2, prod(nr_splines)]));
    assert(Utils.iscloseall(sum(basis_matrix, 2), ones(n_data^2, 1)));
        
    % test if Bspline.fit_Pspline() with lam == 0 produces the same as
    % Bspline.fit()
    [coef_tpbs, basis_matrix_tpbs, knots_tpbs] = Bspline.fit(X,y, nr_splines, orders, knot_types);
    [coef_tpps0, basis_matrix_tpps0, knots_tpps0] = Bspline.fit_Pspline(X,y,0,nr_splines, orders, knot_types);
    assert(Utils.iscloseall(coef_tpbs, coef_tpps0, 1e-8));
    
end

function test_predict_Bspline(testCase)
    rng(2);
    n_data = 250;
    n_data_test = 10;
    x = rand(n_data,1);
    x_test = rand(n_data_test, 1);
    y = sin(6*x)+ 0.2*randn(n_data,1);
    
    nr_splines = 25;
    order = 3;
    knot_type = "e";
    
    [coef, basis_matrix, knots] = Bspline.fit(x,y,nr_splines,order,knot_type);
    
    s = Bspline.predict(x_test, knots, coef, order);
    
    assert(length(s) == n_data_test);
    
end

function test_predict_tensorproduct_Bspline(testCase)
    rng(2); % set seed for random number generator
    n_data = 100;
    n_data_test = 10;
    x1 = rand(100,1);
    x2 = rand(100,1);
    x_test = [rand(n_data_test,1), rand(n_data_test,1)];
    [x1g, x2g] = meshgrid(x1,x2);
    X = [x1g(:), x2g(:)];
    y = sin(X(:,1)*6) + 2*X(:,2).^2 + randn(n_data^2,1)*0.2;
    
    nr_splines = [12,8];
    orders = [3,3];
    knot_types = ["e", "e"];
     
    [coef, basis_matrix, knots] = Bspline.fit(X,y, nr_splines, orders, knot_types);
    
    s = Bspline.predict(x_test, knots, coef, orders);
    
    assert(length(s) == n_data_test);
    
end

function test_calc_GCV(testCase)
    rng(2);
    n_data = 150;
    x = rand(n_data,1);
    y = sin(6*x) + randn(n_data,1)*0.2;
    
    nr_splines = 50;
    order = 3;
    knot_type = "e";
    plot_ = 1;
    nlam = 100;
    
    best_lam = Bspline.calc_GCV(x,y,nlam,plot_, nr_splines, order, knot_type);
    
    [coef, basis_matrix, knots] = Bspline.fit_Pspline(x,y,best_lam, nr_splines);
        
    assert(best_lam > 0);
    assert(all(size(best_lam) == [1,1]));
    
end


function test_calc_GCV_2d(testCase)
    rng(2);
    n_data = 150;
    x1 = rand(n_data,1);
    x2 = rand(n_data,1);
    X = [x1, x2];
    y = sin(6*X(:,1)) + 2*X(:,2).^2 + randn(n_data,1)*0.2;
    
    nr_splines = [12,8];
    order = [3,3];
    knot_type = ["e", "e"];
    plot_ = 1;
    nlam = 100;
    
    best_lam = Bspline.calc_GCV_2d(X,y,nlam,plot_, nr_splines, order, knot_type);
    [coef, basis_matrix, knots] = Bspline.fit_Pspline(X,y,best_lam, nr_splines);
        
    assert(best_lam > 0);
    assert(all(size(best_lam) == [1,1]));
    
end



