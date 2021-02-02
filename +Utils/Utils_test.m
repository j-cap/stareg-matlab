function tests = Utils_test
%%
% Main Function 
%
% run
% >>runtests('Utils_test.m')
tests = functiontests(localfunctions);
end


function test_mapping_matrix_order_1(testCase)
    
    constraint = ["inc", "dec", "peak", "valley"];
    nr_splines = 50;
    for i=1:numel(constraint)
        D = Utils.mapping_matrix(constraint(i), nr_splines);     
        assert(D(1,1) == -1);
        assert(D(1,2) == 1);
        assert(D(end,end-1) == -1);
        assert(D(end,end) == 1);
        assert(Utils.iscloseall(sum(D,2), zeros(nr_splines-1,1)));
    end
end

function test_mapping_matrix_order_2(testCase)
    
    constraint = ["conc", "conv"];
    nr_splines = 50;
    for i=1:numel(constraint)
        D = Utils.mapping_matrix(constraint(i), nr_splines);
        assert(D(1,1) == 1);
        assert(D(1,2) == -2);
        assert(D(1,3) == 1);
        assert(D(end,end-2) == 1);
        assert(D(end,end-1) == -2);
        assert(D(end,end) == 1); 
        assert(Utils.iscloseall(sum(D,2), zeros(nr_splines-2,1)));
    end
end


function test_mapping_matrix_tp_order_1(testCase)
    
    constraint = ["inc", "dec", "peak", "valley"];
    nr_splines = [12,8];
    dim = 1;
    for i=1:numel(constraint)
        D = Utils.mapping_matrix_tp(constraint(i), nr_splines, dim);     
        assert(D(1,1) == -1);
        assert(D(1,2) == 1);
        assert(D(end,end-1) == -1);
        assert(D(end,end) == 1);
        assert(Utils.iscloseall(sum(D,2), zeros(prod(nr_splines)-nr_splines(dim+1),1)));
    end
    
    dim = 2;
    for i=1:numel(constraint)
        D = Utils.mapping_matrix_tp(constraint(i), nr_splines, dim);     
        assert(D(1,1) == -1);
        assert(D(1,1+nr_splines(dim-1)) == 1);
        assert(D(end,end-nr_splines(dim-1)) == -1);
        assert(D(end,end) == 1);
        assert(Utils.iscloseall(sum(D,2), zeros(prod(nr_splines)-nr_splines(dim-1),1)));
    end
     
end

function test_mapping_matrix_tp_order_2(testCase)
    
    constraint = ["conv", "conc"];
    nr_splines = [12,8];
    
    dim = 1; % test for constraint in dimension 1
    for i=1:numel(constraint)
        D = Utils.mapping_matrix_tp(constraint(i), nr_splines, dim);     
        assert(D(1,1) == 1);
        assert(D(1,2) == -2);
        assert(D(1,3) == 1);
        assert(D(end,end-2) == 1);
        assert(D(end,end-1) == -2);
        assert(D(end,end) == 1);
        assert(Utils.iscloseall(sum(D,2), zeros(nr_splines(2) * (nr_splines(dim)-2),1)));
    end
    
    dim = 2; % test for constraint in dimension 2
    for i=1:numel(constraint)
        D = Utils.mapping_matrix_tp(constraint(i), nr_splines, dim);     
        assert(D(1,1) == 1);
        assert(D(1,1+nr_splines(1)) == -2);
        assert(D(1,1+2*nr_splines(1)) == 1);
        assert(D(end,end-2*nr_splines(1)) == 1);
        assert(D(end,end-nr_splines(1)) == -2);
        assert(D(end,end) == 1);
        assert(Utils.iscloseall(sum(D,2), zeros(nr_splines(1) * (nr_splines(dim)-2),1)));
    end
     
    % To plot the matrix, use the commands below
    %>>figure()
    %>>image(D, 'CDataMapping','scaled'); colorbar;
end

function test_check_constraint_inc(testCase)

    x = linspace(0,10,10);
    y = 0;
    B = 0;
    
    constraint = "inc";
    
    coef = x.^2; % test for increasing coefficients -> all zero
    v = Utils.check_constraint(coef, constraint, y, B);
    assert(all(v == zeros(length(x)-1, 1)));
    
    coef = -x.^2; % test for decreasing coefficients -> all one
    v = Utils.check_constraint(coef, constraint, y, B);
    assert(all(v == ones(length(x)-1, 1)));
    
end

function test_check_constraint_dec(testCase)

    x = linspace(0,10,10);
    y = 0;
    B = 0;
    
    constraint = "dec";
    
    coef = -x.^2; % test for decreasing coefficients -> all zero
    v = Utils.check_constraint(coef, constraint, y, B);
    assert(all(v == zeros(length(x)-1, 1)));
    
    coef = x.^2; % test for increasing coefficients -> all one
    v = Utils.check_constraint(coef, constraint, y, B);
    assert(all(v == ones(length(x)-1, 1)));
    
end

function test_check_constraint_conv(testCase)

    x = linspace(0,10,10);
    y = 0;
    B = 0;
    
    constraint = "conv";
    
    coef = exp(x); % test for convex coefficients -> all zero
    v = Utils.check_constraint(coef, constraint, y, B);
    assert(all(v == zeros(length(x)-2, 1)));
    
    coef = x; % test for non-convex coefficients -> all one
    v = Utils.check_constraint(coef, constraint, y, B);
    assert(all(v == ones(length(x)-2, 1)));
    
end

function test_check_constraint_conc(testCase)

    x = linspace(0,10,10);
    y = 0;
    B = 0;
    
    constraint = "conc";
    
    coef = -exp(x); % test for concave coefficients -> all zero
    v = Utils.check_constraint(coef, constraint, y, B);
    assert(all(v == zeros(length(x)-2, 1)));
    
    coef = x; % test for increasing coefficients -> all zero
    v = Utils.check_constraint(coef, constraint, y, B);
    assert(all(v == ones(length(x)-2, 1)));
    
end

function test_check_constraint_peak(testCase)
    
    x = linspace(0,1,100)';
    y = exp(-(x-0.4).^2);
  
    nr_splines = 50;
    order = 3;
    knot_type = "e";
    
    constraint = "peak";

    % test for peak behavior coefficients
    best_lam = Bspline.calc_GCV(x,y,100,0,nr_splines,order,knot_type);
    [coef, basis_matrix, knots] = Bspline.fit_Pspline(x,y,best_lam,nr_splines,order, knot_type);
    v = Utils.check_constraint(coef, constraint, y, basis_matrix);
    assert(all(v == zeros(nr_splines-1,1)));
    
    % TODO: make the peak and valley constraint more robust

end

function test_check_constraint_valley(testCase)
    
    x = linspace(0,1,100)';
    y = -exp(-(x-0.4).^2);
    
    nr_splines = 50;
    order = 3;
    knot_type = "e";
    
    constraint = "valley";

    % test for peak behavior coefficients
    best_lam = Bspline.calc_GCV(x,y,100,0,nr_splines,order,knot_type);
    [coef, basis_matrix, knots] = Bspline.fit_Pspline(x,y,best_lam,nr_splines,order, knot_type);
    v = Utils.check_constraint(coef, constraint, y, basis_matrix);
    assert(all(v == zeros(nr_splines-1,1)));
    
    % TODO: make the peak and valley constraint more robust

end

function test_check_constraint_dim1(testcase)

    nr_splines = [8,12];

    constraint = "inc";
    % test with increasing coef -> v should be all zeros
    coef = linspace(0,1, prod(nr_splines))';
    v = Utils.check_constraint_dim1(coef, constraint, nr_splines);
    assert(all(v == zeros((nr_splines(1)-1)*nr_splines(2),1)))
    
    % test with decreasing coef -> v should be all ones
    coef = -linspace(0,1, prod(nr_splines))';
    v = Utils.check_constraint_dim1(coef, constraint, nr_splines);
    assert(all(v == ones((nr_splines(1)-1)*nr_splines(2),1)))
    
    constraint = "dec";
    % test with decreasing coef -> v should be all zeros
    coef = -linspace(0,1, prod(nr_splines))';
    v = Utils.check_constraint_dim1(coef, constraint, nr_splines);
    assert(all(v == zeros((nr_splines(1)-1)*nr_splines(2),1)))
    
    % test with increasing coef -> v should be all ones
    coef = linspace(0,1, prod(nr_splines))';
    v = Utils.check_constraint_dim1(coef, constraint, nr_splines);
    assert(all(v == ones((nr_splines(1)-1)*nr_splines(2),1)))
    
end


function test_check_constraint_dim2(testcase)

    nr_splines = [8,12];

    constraint = "inc";
    % test with increasing coef -> v should be all zeros
    coef = linspace(0,1, prod(nr_splines))';
    v = Utils.check_constraint_dim2(coef, constraint, nr_splines);
    assert(all(v == zeros((nr_splines(2)-1)*nr_splines(1),1)))
    
    % test with decreasing coef -> v should be all ones
    coef = -linspace(0,1, prod(nr_splines))';
    v = Utils.check_constraint_dim2(coef, constraint, nr_splines);
    assert(all(v == ones((nr_splines(2)-1)*nr_splines(1),1)))
    
    constraint = "dec";
    % test with decreasing coef -> v should be all zeros
    coef = -linspace(0,1, prod(nr_splines))';
    v = Utils.check_constraint_dim2(coef, constraint, nr_splines);
    assert(all(v == zeros((nr_splines(2)-1)*nr_splines(1),1)))
    
    % test with increasing coef -> v should be all ones
    coef = linspace(0,1, prod(nr_splines))';
    v = Utils.check_constraint_dim2(coef, constraint, nr_splines);
    assert(all(v == ones((nr_splines(2)-1)*nr_splines(1),1)))
    
end

function test_check_constraint_full(testCase)
    rng(2);
    n_data = 500;
    X = [rand(n_data,1), rand(n_data,2)];
    y = 4*X(:,1) + X(:,1) .* X(:,2) + randn(n_data,1)*0.2;

    description = {
        ["s(1)", "inc", 50, 3000, "e"]; 
        ["t(1,2)", "none,inc", "12,20", "2000,2000", "e,q"]};
    model = Stareg.create_model_from_description(description, X, y);  
    [B, S, K, W, coef] = Stareg.create_model_matrices(model);

    [w, W_compare] = Utils.check_constraint_full(coef, model, B, y);
    
    assert(length(w.f1.v) == 50-1);
    assert(length(w.f2.v1) == 12*20);
    assert(length(w.f2.v2) == 12*19);
    assert(length(W_compare) == 50-1 + 12*20 + 12*19); 
    
end

function test_create_constraint_matrix(testCase)

    rng(2);
    n_data = 500;
    X = [rand(n_data,1), rand(n_data,2)];
    y = 4*X(:,1) + X(:,1) .* X(:,2) + randn(n_data,1)*0.2;

    description = {
        ["s(1)", "inc", 50, 3000, "e"]; 
        ["t(1,2)", "none,inc", "12,20", "2000,2000", "e,q"]};
    model = Stareg.create_model_from_description(description, X, y);  
    [B, S, K, W, coef] = Stareg.create_model_matrices(model);

    [w, W_compare] = Utils.check_constraint_full(coef, model, B, y);
    
    KK = Utils.create_constraint_matrix(model, w);
    
    assert(all(size(KK) == [50+12*20, 50+12*20]));
    
    
end

