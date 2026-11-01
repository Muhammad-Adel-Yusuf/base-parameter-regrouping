%% ========================================================================
% B_Creating_A_matrix_base.m  (ABB IRB 120)
%
% Automatic base-parameter regrouping (Algorithm 1 of the paper) with the
% link dimensions kept symbolic. See regroup_base_parameters.m.
%
% INPUT : torque_equations/Torque_1.txt ... Torque_6.txt
% OUTPUT: results/base_regressor.mat, results/theta_base.txt,
%         results/theta_base_list.txt, results/A_matrix_base.txt,
%         results/regrouping_summary.txt
%% ========================================================================
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here);

opt.robot_name      = 'IRB120';
opt.torque_folder   = fullfile(here, 'torque_equations');
opt.results_folder  = fullfile(here, 'results');
% link dimensions (symbolic in the model; numeric values are used for the
% classical-count comparison and for B2_known_geometry_regrouping.m)
opt.geometry_names  = {'L1','L2','L3','L4','L5','L6'};
opt.geometry_values = [0.29 0.27 0.07 0.302 0.072 0];
opt.known_geometry  = false;

summary = regroup_base_parameters(opt);
