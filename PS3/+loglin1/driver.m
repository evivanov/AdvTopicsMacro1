%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'loglin1';
M_.dynare_version = '4.6.4';
oo_.dynare_version = '4.6.4';
options_.dynare_version = '4.6.4';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('loglin1.log');
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'eps'};
M_.exo_names_tex(1) = {'eps'};
M_.exo_names_long(1) = {'eps'};
M_.endo_names = cell(2,1);
M_.endo_names_tex = cell(2,1);
M_.endo_names_long = cell(2,1);
M_.endo_names(1) = {'k'};
M_.endo_names_tex(1) = {'k'};
M_.endo_names_long(1) = {'k'};
M_.endo_names(2) = {'z'};
M_.endo_names_tex(2) = {'z'};
M_.endo_names_long(2) = {'z'};
M_.endo_partitions = struct();
M_.param_names = cell(9,1);
M_.param_names_tex = cell(9,1);
M_.param_names_long = cell(9,1);
M_.param_names(1) = {'alpha'};
M_.param_names_tex(1) = {'alpha'};
M_.param_names_long(1) = {'alpha'};
M_.param_names(2) = {'beta'};
M_.param_names_tex(2) = {'beta'};
M_.param_names_long(2) = {'beta'};
M_.param_names(3) = {'rho'};
M_.param_names_tex(3) = {'rho'};
M_.param_names_long(3) = {'rho'};
M_.param_names(4) = {'sigma'};
M_.param_names_tex(4) = {'sigma'};
M_.param_names_long(4) = {'sigma'};
M_.param_names(5) = {'A'};
M_.param_names_tex(5) = {'A'};
M_.param_names_long(5) = {'A'};
M_.param_names(6) = {'B'};
M_.param_names_tex(6) = {'B'};
M_.param_names_long(6) = {'B'};
M_.param_names(7) = {'C'};
M_.param_names_tex(7) = {'C'};
M_.param_names_long(7) = {'C'};
M_.param_names(8) = {'D'};
M_.param_names_tex(8) = {'D'};
M_.param_names_long(8) = {'D'};
M_.param_names(9) = {'E'};
M_.param_names_tex(9) = {'E'};
M_.param_names_long(9) = {'E'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 2;
M_.param_nbr = 9;
M_.orig_endo_nbr = 2;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
options_.linear = true;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
options_.linear_decomposition = false;
M_.nonzero_hessian_eqs = [];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.orig_eq_nbr = 2;
M_.eq_nbr = 2;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 1 3 5;
 2 4 6;]';
M_.nstatic = 0;
M_.nfwrd   = 0;
M_.npred   = 0;
M_.nboth   = 2;
M_.nsfwrd   = 2;
M_.nspred   = 2;
M_.ndynamic   = 2;
M_.dynamic_tmp_nbr = [0; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , '1' ;
  2 , 'name' , 'z' ;
};
M_.mapping.k.eqidx = [1 ];
M_.mapping.z.eqidx = [1 2 ];
M_.mapping.eps.eqidx = [2 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [1 2 ];
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(2, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(9, 1);
M_.endo_trends = struct('deflator', cell(2, 1), 'log_deflator', cell(2, 1), 'growth_factor', cell(2, 1), 'log_growth_factor', cell(2, 1));
M_.NNZDerivatives = [8; 0; -1; ];
M_.static_tmp_nbr = [0; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(2) = 0.99;
beta = M_.params(2);
M_.params(1) = 0.33;
alpha = M_.params(1);
M_.params(3) = 0.95;
rho = M_.params(3);
M_.params(4) = 0.01;
sigma = M_.params(4);
M_.params(5) = M_.params(1)*M_.params(2);
A = M_.params(5);
M_.params(6) = M_.params(2)*(-(M_.params(1)^2))-1;
B = M_.params(6);
M_.params(7) = M_.params(1);
C = M_.params(7);
M_.params(8) = M_.params(2)*(-M_.params(1));
D = M_.params(8);
M_.params(9) = 1;
E = M_.params(9);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
model_diagnostics(M_,options_,oo_);
options_.order = 1;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
save('loglin1_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('loglin1_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('loglin1_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('loglin1_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('loglin1_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('loglin1_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('loglin1_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
