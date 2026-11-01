%% ========================================================================
% B4_actuation_known_geometry.m  (ABB IRB 120)
%
% Automatic base-parameter regrouping with the joint-side actuation terms
% added to every joint torque as symbolic parameters:
%     tau_j <- tau_j + Ia_j*qdd_j + Fv_j*qd_j + Fc_j*sign(qd_j)
% (rotor inertia reflected to the joint, viscous and Coulomb friction).
% Link dimensions substituted (as in B2_known_geometry_regrouping.m).
% No numerical values of these terms are assumed: the regrouping and its
% verification use random values for every symbol.
% Results go to results/actuation_known_geometry/.
%% ========================================================================
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here);

opt.robot_name      = 'IRB120';
opt.torque_folder   = fullfile(here, 'torque_equations');
opt.results_folder  = fullfile(here, 'results', 'actuation_known_geometry');
opt.geometry_names  = {'L1','L2','L3','L4','L5','L6'};
opt.geometry_values = [0.29 0.27 0.07 0.302 0.072 0];
opt.known_geometry  = true;
opt.actuation       = true;

summary = regroup_base_parameters(opt);
