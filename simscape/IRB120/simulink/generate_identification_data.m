%% ========================================================================
% generate_identification_data.m  (ABB IRB 120)
%
% 1) Builds the Simscape Multibody model IRB120_identification.slx from the
%    manufacturer description shipped with MATLAB (loadrobot abbIrb120:
%    link geometry, masses, centres of mass and inertia tensors).
% 2) Drives every joint with the quintic-ramped multisine excitation and
%    lets Simscape compute the joint torques.
% 3) Saves ../data/simscape_data.mat in the MDH joint convention.
% 4) Writes ../nominal_parameters_IRB120.m with the link data expressed in
%    the MDH frames of A_Full_Dynamics_IRB120_parallel_axis.m.
%
% Requires: Simulink, Simscape Multibody, Robotics System Toolbox.
%% ========================================================================

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(here);

robot_name = 'abbIrb120';
bodies = {'link_1','link_2','link_3','link_4','link_5','link_6'};

% Joint convention of the symbolic model coincides with the manufacturer's
S     = ones(1,6);
delta = zeros(1,6);

% Link dimensions of the MDH table (A_Full_Dynamics_IRB120_parallel_axis.m)
geom = struct('L1', 0.29, 'L2', 0.27, 'L3', 0.07, 'L4', 0.302, 'L5', 0.072, 'L6', 0);
mdh.alpha  = [0 -pi/2 0 -pi/2 pi/2 -pi/2];
mdh.a      = [0 0 geom.L2 geom.L3 0 0];
mdh.d      = [geom.L1 0 0 geom.L4 0 geom.L5];
mdh.theta0 = [0 -pi/2 0 0 0 0];

r = loadrobot(robot_name, 'DataFormat', 'row');
lim = cell2mat(cellfun(@(b) b.Joint.PositionLimits, r.Bodies(2:7).', 'UniformOutput', false));

cfg.robot      = robot_name;
cfg.model_name = 'IRB120_identification';
cfg.label      = 'IRB 120';
cfg.folder     = here;
cfg.S = S;  cfg.delta = delta;
cfg.q_min      = lim(:,1).';
cfg.q_max      = lim(:,2).';
cfg.v_max      = deg2rad([250 250 250 320 320 420]);       % IRB 120 axis speeds
cfg.frac       = 0.8;
cfg.f          = [0.050 0.065 0.085 0.105 0.125 0.145];
cfg.phi        = [0 pi/3 pi/5 pi/7 pi/4 pi/6];
cfg.T_ramp     = 5;
cfg.T_end      = 40;
cfg.Ts         = 0.004;

data = build_and_simulate(cfg);
save(fullfile(root, 'data', 'simscape_data.mat'), '-struct', 'data');
fprintf('Saved %d samples to data/simscape_data.mat\n', numel(data.time_data));
fprintf('Limits (deg): %s\n', mat2str(rad2deg(lim), 4));
fprintf('Amplitudes (deg): %s\n', mat2str(rad2deg(data.amplitude), 4));
fprintf('Max |qd| / v_max: %s\n', mat2str(max(abs(data.manufacturer.qd))./cfg.v_max, 3));
fprintf('Max |tau| (N m) : %s\n', mat2str(max(abs(data.T_data)), 3));

P = extract_nominal_parameters(robot_name, bodies, mdh, S, delta);
fprintf('Body-to-MDH attachment error: %.2e\n', P.attach_error);
write_nominal_file(fullfile(root, 'nominal_parameters_IRB120.m'), 'IRB120', robot_name, P, geom);
