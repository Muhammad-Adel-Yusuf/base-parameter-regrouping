%% ========================================================================
% B2_known_geometry_regrouping.m  (ABB IRB 120)
%
% Same regrouping as B_Creating_A_matrix_base.m, but with the link
% dimensions substituted before extraction (classical setting, geometry
% known from kinematic calibration). Results go to results/known_geometry/.
%% ========================================================================
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here);

opt.robot_name      = 'IRB120';
opt.torque_folder   = fullfile(here, 'torque_equations');
opt.results_folder  = fullfile(here, 'results', 'known_geometry');
opt.geometry_names  = {'L1','L2','L3','L4','L5','L6'};
opt.geometry_values = [0.29 0.27 0.07 0.302 0.072 0];
opt.known_geometry  = true;

summary = regroup_base_parameters(opt);
