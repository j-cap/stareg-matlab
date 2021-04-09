
%% Script to run all test at once

% run Bspline package test
res_Bspline = run(Bspline.Bspline_test);
% run Utils package test
res_Utils = run(Utils.Utils_test);
% run Stareg package test
res_Stareg = run(Stareg.Stareg_test);

res_Bspline
res_Utils
res_Stareg


