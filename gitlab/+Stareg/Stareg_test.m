function tests = Stareg_test
%%
% Main Function 
tests = functiontests(localfunctions);
end

function test_predict(testCase)
    rng(2);
    description = {
        ["s(1)", 100, "inc", 3000, "e"]; 
        ["s(2)", 10, "peak", 300, "e"];
        ["t(1,2)", "12,20", "inc,none", "2000,2000", "e,q"]
        };
    n_data = 500;
    X = [rand(n_data,1), rand(n_data,1)];
    y = 4*X(:,1) + exp(-(X(:,2)-0.4).^2 ./ 0.01) + X(:,1) .* X(:,2) + randn(n_data,1)*0.2;
    
    [coef, basis_matrix, model] = Stareg.fit(description, X, y);
    n_pred = 100;
    Xpred = [linspace(0,1,n_pred)', rand(n_pred, 1)];
    s = Stareg.predict(Xpred, model, coef);
    
    assert(length(s) == n_pred); 

end



function test_create_model_from_description(testCase)
    rng(2);
    description = {
        ["s(1)", 100, "inc", 3000, "e"]; 
        ["s(2)", 10, "peak", 300, "e"];
        ["t(1,2)", "12,20", "inc,none", "2000,2000", "e,q"]
        };
    n_data = 500;
    X = [rand(n_data,1), rand(n_data,2)];
    y = 4*X(:,1) + exp(-(X(:,2)-0.4).^2 ./ 0.01) + X(:,1) .* X(:,2) + randn(n_data,1)*0.2;

    model = Stareg.create_model_from_description(description, X, y);
    
    assert(length(fieldnames(model)) == 3);
    assert(model.f1.type == "s(1)");
    assert(model.f1.nr_splines == 100);
    assert(model.f1.constraint == "inc");
    assert(model.f1.lam_c == 3000);
    assert(model.f1.knot_type == "e");
    assert(all(size(model.f1.B) == [n_data, 100]));
    assert(all(size(model.f1.Ds) == [100-2, 100]));
    assert(all(size(model.f1.Dc) == [100-1, 100]));
    assert(strcmp(class(model.f1.best_lam), class(1.0)));
    assert(length(model.f1.coef_pls) == 100);
    assert(length(model.f1.v) == 100-1);
    
    assert(model.f2.type == "s(2)");
    assert(model.f2.nr_splines == 10);
    assert(model.f2.constraint == "peak");
    assert(model.f2.lam_c == 300);
    assert(model.f2.knot_type == "e");
    assert(all(size(model.f2.B) == [n_data, 10]));
    assert(all(size(model.f2.Ds) == [10-2, 10]));
    assert(all(size(model.f2.Dc) == [10-1, 10]));
    assert(strcmp(class(model.f2.best_lam), class(1.0)));
    assert(length(model.f2.coef_pls) == 10);
    assert(length(model.f2.v) == 10-1);
    
    assert(model.f3.type == "t(1,2)");
    assert(all(diag(model.f3.nr_splines == [12,20])));
    assert(all(diag(model.f3.constraint == ["inc", "none"])));
    assert(all(diag(model.f3.lam_c == [2000,2000])));
    assert(all(diag(model.f3.knot_type == ["e","q"])));
    assert(all(size(model.f3.B) == [n_data, 240]));
    assert(all(size(model.f3.Ds.Ds1) == [10*20, 240]));
    assert(all(size(model.f3.Ds.Ds2) == [12*18, 240]));
    assert(all(size(model.f3.Dc.Dc1) == [11*20, 240]));
    assert(all(size(model.f3.Dc.Dc2) == [12*20, 240]));

    assert(strcmp(class(model.f3.best_lam), class(1.0)));
    assert(length(model.f3.coef_pls) == 12*20);
    assert(length(model.f3.v.v1) == 11*20);
    assert(length(model.f3.v.v2) == 12*20);
    
end

function test_create_model_matrices(testCase)
    rng(2);
    description = {
        ["s(1)", 100, "inc", 3000, "e"]; 
        ["s(2)", 10, "conc", 300, "e"];
        ["t(1,2)", "12,20", "inc,none", "2000,2000", "e,q"]
        };
    n_data = 500;
    X = [rand(n_data,1), rand(n_data,2)];
    y = 4*X(:,1) + exp(-(X(:,2)-0.4).^2 ./ 0.01) + X(:,1) .* X(:,2) + randn(n_data,1)*0.2;
        
    model = Stareg.create_model_from_description(description, X, y);
    [B, S, K, W, coef] = Stareg.create_model_matrices(model);
    
    assert(all(size(B) == [n_data, 100+10+12*20]));
    assert(all(size(S) == [100 + 10 + 12*20, 100 + 10 + 12*20]));
    assert(all(size(K) == [100 + 10 + 12*20, 100 + 10 + 12*20]));
    assert(length(W) == 100-1 + 10-2 + (12-1)*20 + 12*20);
    assert(length(coef) == 100+10+12*20);

end

function test_fit_1d_inc(testCase)
    rng(2);
    n_data = 100;
    x = rand(n_data,1);
    y = sin(6*x) + 4*x.^2 + randn(n_data, 1)*0.2;

    description = {["s(1)", 80, "inc", 6000, "e"]; };
    [coef, basis_matrix, model, reduced_model, coef_list] = Stareg.fit(description, x, y);

    fn = fieldnames(coef_list);
    assert(length(coef) == 80);
    assert(length(coef_list.(fn{end})) == 80);
    assert(all(size(basis_matrix) == [100,80]));
    
    % plot the iteration
    %figure();
    %scatter(x,y, 'DisplayName','Data'); hold on;
    %for i=1:numel(fn)
    %    scatter(x, basis_matrix*coef.(fn{i}), 'x', 'DisplayName', string(i))
    %end
    %grid(); legend();
    
end

function test_fit_1d_dec(testCase)
    rng(2);
    n_data = 100;
    x = rand(n_data,1);
    y = sin(6*x) + 4*x.^2 + randn(n_data, 1)*0.2;

    description = {["s(1)", 80, "dec", 6000, "e"]; };
    [coef, basis_matrix, model, reduced_model, coef_list] = Stareg.fit(description, x, y);

    fn = fieldnames(coef_list);
    assert(length(coef) == 80);
    assert(length(coef_list.(fn{end})) == 80);
    assert(all(size(basis_matrix) == [100,80]));
    
    % plot the iteration
    %figure();
    %scatter(x,y, 'DisplayName','Data'); hold on;
    %for i=1:numel(fn)
    %    scatter(x, basis_matrix*coef.(fn{i}), 'x', 'DisplayName', string(i))
    %end
    %grid(); legend();
    
end

function test_fit_1d_none(testCase)
    rng(2);
    n_data = 100;
    x = rand(n_data,1);
    y = sin(6*x) + 4*x.^2 + randn(n_data, 1)*0.2;

    description = {["s(1)", 80, "none", 6000, "e"]; };
    [coef, basis_matrix, model, reduced_model, coef_list] = Stareg.fit(description, x, y);

    fn = fieldnames(coef_list);
    assert(length(coef) == 80);
    assert(length(coef_list.(fn{end})) == 80);
    assert(all(size(basis_matrix) == [100,80]));
    
    % plot the iteration
    %figure();
    %scatter(x,y, 'DisplayName','Data'); hold on;
    %for i=1:numel(fn)
    %    scatter(x, basis_matrix*coef.(fn{i}), 'x', 'DisplayName', string(i))
    %end
    %grid(); legend();
    
end

function test_fit_1d_peak(testCase)
    rng(2);
    n_data = 100;
    x = rand(n_data,1);
    y = 4*exp(-(x-0.4).^2 / 0.01) + x.^2 + randn(n_data, 1)*0.2;
    
    description = {["s(1)", 80, "peak", 6000, "e"]; };
    [coef, basis_matrix, model, reduced_model, coef_list] = Stareg.fit(description, x, y);

    fn = fieldnames(coef_list);
    assert(length(coef) == 80);
    assert(length(coef_list.(fn{end})) == 80);
    assert(all(size(basis_matrix) == [100,80]));
    
    % plot the iteration
    %figure();
    %scatter(x,y, 'DisplayName','Data'); hold on;
    %for i=1:numel(fn)
    %    scatter(x, basis_matrix*coef.(fn{i}), 'x', 'DisplayName', string(i))
    %end
    %grid(); legend();
    
end

function test_fit_1d_valley(testCase)
    rng(2);
    n_data = 100;
    x = rand(n_data,1);
    y = -4*exp(-(x-0.4).^2 / 0.01) + x.^2 + randn(n_data, 1)*0.2;
    
    description = {["s(1)", 80, "valley", 6000, "e"]; };
    [coef, basis_matrix, model, reduced_model, coef_list] = Stareg.fit(description, x, y);

    fn = fieldnames(coef_list);
    assert(length(coef) == 80);
    assert(length(coef_list.(fn{end})) == 80);
    assert(all(size(basis_matrix) == [100,80]));
    
    % plot the iteration
    %figure();
    %scatter(x,y, 'DisplayName','Data'); hold on;
    %for i=1:numel(fn)
    %    scatter(x, basis_matrix*coef.(fn{i}), 'x', 'DisplayName', string(i))
    %end
    %grid(); legend();
    
end


function test_fit_2d(testCase)
    rng(2);
    n_data = 5000;
    X = [rand(n_data,1), rand(n_data,1)];
    y = 4*exp(-(X(:,1)-0.4).^2 / 0.01) + X(:,2).^2 + randn(n_data, 1)*0.2;
    
    description = {
        ["s(1)", 80, "peak", 3000, "e"]; 
        ["t(1,2)", "12,8", "none,inc", "2000,2000", "e,q"]
        };
    
    [coef, basis_matrix, model, reduced_model, coef_list] = Stareg.fit(description, X, y);

    fn = fieldnames(coef_list);
    assert(length(coef) == 80+12*8);
    assert(length(coef_list.(fn{end})) == 80+12*8);
    assert(all(size(basis_matrix) == [n_data,80+12*8]));
    
    % plot the iteration
    %figure();
    %scatter3(X(:,1),X(:,2), y, 'DisplayName','Data'); hold on;
    %for i=1:numel(fn)
    %    scatter3(X(:,1), X(:,2), basis_matrix*coef.(fn{i}), 'x', 'DisplayName', string(i))
    %end
    %grid(); legend();
        
end




