%% ========================================================================
% generate_identification_data.m  (Franka Emika Panda)
%
% 1) Builds the Simscape Multibody model Panda_identification.slx from
%    panda_robot.m: link geometry and meshes of loadrobot frankaEmikaPanda,
%    masses, centres of mass and inertia tensors of Gaz et al. (2019).
% 2) Drives every joint with the quintic-ramped multisine excitation and
%    lets Simscape compute the joint torques.
% 3) Saves ../data/simscape_data.mat in the MDH joint convention.
% 4) Writes ../nominal_parameters_Panda.m with the link data expressed in
%    the MDH frames of A_Full_Dynamics_Panda_parallel_axis.m.
%
% Requires: Simulink, Simscape Multibody, Robotics System Toolbox.
%% ========================================================================

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(here);  addpath(root);

robot_name = 'frankaEmikaPanda';
bodies = arrayfun(@(k) sprintf('panda_link%d', k), 1:7, 'UniformOutput', false);

% Joint convention of the symbolic model coincides with the manufacturer's
S     = ones(1,7);
delta = zeros(1,7);

% Link dimensions of the MDH table (A_Full_Dynamics_Panda_parallel_axis.m)
geom = struct('L1', 0.333, 'L2', 0.316, 'L3', 0.0825, 'L4', 0.0825, 'L5', 0.384, 'L6', 0.088);
mdh.alpha  = [0 -pi/2 pi/2 pi/2 -pi/2 pi/2 pi/2];
mdh.a      = [0 0 0 geom.L3 -geom.L4 0 geom.L6];
mdh.d      = [geom.L1 0 geom.L2 0 geom.L5 0 0];
mdh.theta0 = zeros(1,7);

r = panda_robot();
lim = cell2mat(cellfun(@(b) b.Joint.PositionLimits, r.Bodies(1:7).', 'UniformOutput', false));

cfg.robot      = robot_name;
cfg.model_name = 'Panda_identification';
cfg.label      = 'Panda';
cfg.folder     = here;
cfg.S = S;  cfg.delta = delta;
cfg.q_min      = lim(:,1).';
cfg.q_max      = lim(:,2).';
cfg.v_max      = [2.175 2.175 2.175 2.175 2.61 2.61 2.61]; % Panda joint speed limits (rad/s)
cfg.frac       = 0.8;
cfg.f          = [0.050 0.065 0.085 0.105 0.125 0.145 0.165];
cfg.phi        = [0 pi/3 pi/5 pi/7 pi/4 pi/6 pi/8];
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
write_nominal_file(fullfile(root, 'nominal_parameters_Panda.m'), 'Panda', robot_name, P, geom);
