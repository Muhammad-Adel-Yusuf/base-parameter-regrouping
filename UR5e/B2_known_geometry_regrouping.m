%% ========================================================================
% B2_known_geometry_regrouping.m  (UR5e)
%
% Same regrouping as B_Creating_A_matrix_base.m, but with the link
% dimensions substituted before extraction (classical setting, geometry
% known from kinematic calibration). Results go to results/known_geometry/.
%% ========================================================================
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here);

opt.robot_name      = 'UR5e';
opt.torque_folder   = fullfile(here, 'torque_equations');
opt.results_folder  = fullfile(here, 'results', 'known_geometry');
opt.geometry_names  = {'L1','L2','L3','L4','L5','L6','LL1','LL2'};
opt.geometry_values = [0.1625 0.425 0.3922 -0.1333 -0.0997 -0.0996 0 0];
opt.known_geometry  = true;

summary = regroup_base_parameters(opt);
