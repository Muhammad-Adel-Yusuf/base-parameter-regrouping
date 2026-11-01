%% ========================================================================
% B3_actuation_regrouping.m  (Franka Emika Panda)
%
% Automatic base-parameter regrouping with the joint-side actuation terms
% added to every joint torque as symbolic parameters:
%     tau_j <- tau_j + Ia_j*qdd_j + Fv_j*qd_j + Fc_j*sign(qd_j)
% (rotor inertia reflected to the joint, viscous and Coulomb friction).
% Link dimensions symbolic (as in B_Creating_A_matrix_base.m).
% No numerical values of these terms are assumed: the regrouping and its
% verification use random values for every symbol.
% Results go to results/actuation/.
%% ========================================================================
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here);

opt.robot_name      = 'Panda';
opt.nDOF            = 7;
opt.torque_folder   = fullfile(here, 'torque_equations');
opt.results_folder  = fullfile(here, 'results', 'actuation');
opt.geometry_names  = {'L1','L2','L3','L4','L5','L6'};
opt.geometry_values = [0.333 0.316 0.0825 0.0825 0.384 0.088];
opt.known_geometry  = false;
opt.actuation       = true;

% Stage 1 on the nested Newton-Euler expressions (A2_regressor_columns_Panda.m)
opt.column_file     = fullfile(here, 'torque_equations', 'regressor_columns.mat');
opt.torque_function = @panda_rnea_numeric;

summary = regroup_base_parameters(opt);
