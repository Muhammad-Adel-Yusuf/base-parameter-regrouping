%% ========================================================================
% generate_identification_data.m  (UR5e)
%
% 1) Builds the Simscape Multibody model UR5e_identification.slx from the
%    manufacturer description shipped with MATLAB (loadrobot universalUR5e:
%    link geometry, masses, centres of mass and inertia tensors).
% 2) Drives every joint with the quintic-ramped multisine excitation and
%    lets Simscape compute the joint torques.
% 3) Saves ../data/simscape_data.mat in the MDH joint convention.
% 4) Writes ../nominal_parameters_UR5e.m with the link data expressed in the
%    MDH frames of A_Full_Dynamics_UR5e_parallel_axis.m.
%
% Requires: Simulink, Simscape Multibody, Robotics System Toolbox.
% Run from a working folder with a short path (Windows 260-character limit
% for Simulink build files).
%% ========================================================================

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(here);

robot_name = 'universalUR5e';
bodies = {'shoulder_link','upper_arm_link','forearm_link','wrist_1_link','wrist_2_link','wrist_3_link'};

% Joint convention of the symbolic model: q = S*q_m + delta, tau = S*tau_m
S     = [1 1 1 -1 -1 -1];
delta = [0 pi/2 0 pi/2 0 0];

% Link dimensions of the MDH table (A_Full_Dynamics_UR5e_parallel_axis.m)
geom = struct('L1', 0.1625, 'L2', 0.425, 'L3', 0.3922, 'L4', -0.1333, ...
              'L5', -0.0997, 'L6', -0.0996, 'LL1', 0, 'LL2', 0);
mdh.alpha  = [0 pi/2 0 -pi -pi/2 pi/2];
mdh.a      = [0 0 geom.L2 geom.L3 0 0];
mdh.d      = [geom.L1 geom.LL1 -geom.LL2 geom.L4 geom.L5 geom.L6];
mdh.theta0 = [-pi pi/2 0 -pi/2 0 pi/2];

% Excitation (manufacturer joint coordinates)
cfg.robot      = robot_name;
cfg.model_name = 'UR5e_identification';
cfg.label      = 'UR5e';
cfg.folder     = here;
cfg.S = S;  cfg.delta = delta;
cfg.q_min      = deg2rad([-170 -120 -170 -180 -120 -360]);
cfg.q_max      = deg2rad([ 170  120  170  180  120  360]);
cfg.v_max      = deg2rad([180 180 180 180 180 180]);       % UR5e joint speed limits
cfg.frac       = 0.8;
cfg.f          = [0.050 0.065 0.085 0.105 0.125 0.145];
cfg.phi        = [0 pi/3 pi/5 pi/7 pi/4 pi/6];
cfg.T_ramp     = 5;
cfg.T_end      = 40;
cfg.Ts         = 0.004;

data = build_and_simulate(cfg);
save(fullfile(root, 'data', 'simscape_data.mat'), '-struct', 'data');
fprintf('Saved %d samples to data/simscape_data.mat\n', numel(data.time_data));
fprintf('Amplitudes (deg): %s\n', mat2str(rad2deg(data.amplitude), 4));
fprintf('Max |qd| / v_max: %s\n', mat2str(max(abs(data.manufacturer.qd))./cfg.v_max, 3));
fprintf('Max |tau| (N m) : %s\n', mat2str(max(abs(data.T_data)), 3));

P = extract_nominal_parameters(robot_name, bodies, mdh, S, delta);
fprintf('Body-to-MDH attachment error: %.2e\n', P.attach_error);
write_nominal_file(fullfile(root, 'nominal_parameters_UR5e.m'), 'UR5e', robot_name, P, geom);
