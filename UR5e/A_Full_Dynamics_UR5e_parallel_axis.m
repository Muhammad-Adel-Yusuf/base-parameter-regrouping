%% ========================================================================
% A_Full_Dynamics_UR5e_parallel_axis.m
%
% Symbolic Euler-Lagrange inverse dynamics of the UR5e (modified DH,
% Khalil-Kleinfinger convention) using ORIGIN-REFERENCED inertia tensors
% (parallel-axis formulation).
%
%   K_i = 1/2 m_i v_i'v_i + m_i v_i'(w_i x c_i) + 1/2 w_i' I_oi w_i
%   V_i = -m_i g0'(R_i c_i + o_i)
%
% where I_oi = I_ci + m_i (c_i'c_i E3 - c_i c_i') is the inertia tensor of
% link i about the ORIGIN of frame {i} (Khalil's XX_i..ZZ_i). In this file
% the symbols Ixx_i..Izz_i therefore denote ORIGIN-referenced inertia.
% The torques are then linear in the 10 standard parameters per link
% [XX XY XZ YY YZ ZZ MX MY MZ M] with MX_i = m_i*Lc_ix, etc.
%
% g is substituted numerically (known constant). Link dimensions L1..L6,
% LL1, LL2 are kept symbolic.
%
% OUTPUT: torque_equations/Torque_1.txt ... Torque_6.txt
% Run time: tens of minutes (symbolic simplification).
%% ========================================================================

clear; clc;

here = fileparts(mfilename('fullpath'));
out_folder = fullfile(here, 'torque_equations');
if ~isfolder(out_folder), mkdir(out_folder); end

syms alpha a theta d
syms th1(t) th2(t) th3(t) th4(t) th5(t) th6(t) real
syms L1 L2 L3 L4 L5 L6 LL1 LL2 real
syms g real

%% Modified DH table (UR5e)  [alpha_{i-1}  a_{i-1}  theta_i  d_i]
%  LL1 and LL2 are the shoulder/elbow lateral offsets (L_e1, L_e2 in the paper)
N = 7;
DH_Table = [ 0      0    th1(t)-pi     L1;
             pi/2   0    th2(t)+pi/2   LL1;
             0      L2   th3(t)       -LL2;
            -pi     L3   th4(t)-pi/2   L4;
            -pi/2   0    th5(t)        L5;
             pi/2   0    th6(t)+pi/2   L6;
             0      0    0             0 ];

T_DH = [cos(theta)              -sin(theta)             0            a;
        sin(theta)*cos(alpha)    cos(theta)*cos(alpha) -sin(alpha) -d*sin(alpha);
        sin(theta)*sin(alpha)    cos(theta)*sin(alpha)  cos(alpha)  d*cos(alpha);
        0                        0                      0            1];

A = cell(1,N);
for i = 1:N
    alpha = DH_Table(i,1); a = DH_Table(i,2); theta = DH_Table(i,3); d = DH_Table(i,4);
    A{i} = subs(T_DH);
end

%% Link symbols
nL = 6;
th = {th1, th2, th3, th4, th5, th6};
m  = sym('m',  [1 nL], 'real');
Lc = sym(zeros(3,nL));
I  = cell(1,nL);
for i = 1:nL
    Lc(:,i) = [sym(sprintf('Lc%dx',i),'real'); sym(sprintf('Lc%dy',i),'real'); sym(sprintf('Lc%dz',i),'real')];
    Ixx = sym(sprintf('Ixx%d',i),'real'); Iyy = sym(sprintf('Iyy%d',i),'real'); Izz = sym(sprintf('Izz%d',i),'real');
    Ixy = sym(sprintf('Ixy%d',i),'real'); Ixz = sym(sprintf('Ixz%d',i),'real'); Iyz = sym(sprintf('Iyz%d',i),'real');
    I{i} = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz];   % ORIGIN-referenced tensor
end

%% Velocity propagation and energies (origin-referenced / parallel-axis form)
w = sym([0;0;0]);  v = sym([0;0;0]);
T0 = sym(eye(4));
KE = sym(0);  PE = sym(0);
for i = 1:nL
    R_im1_i = A{i}(1:3,1:3);  p_im1_i = A{i}(1:3,4);
    R_i_im1 = R_im1_i.';
    th_dot  = diff(th{i}, t);
    v = R_i_im1*(v + cross(w, p_im1_i));        % origin velocity of {i} in {i}
    w = R_i_im1*w + [0;0;th_dot];                % angular velocity of {i} in {i}
    T0 = simplify(T0*A{i});
    R0i = T0(1:3,1:3);  o0i = T0(1:3,4);
    KE = KE + 0.5*m(i)*(v.'*v) + m(i)*(v.'*cross(w, Lc(:,i))) + 0.5*(w.'*I{i}*w);
    PE = PE + m(i)*[0 0 g]*(R0i*Lc(:,i) + o0i);
end
Lag = KE - PE;

%% Euler-Lagrange equations
torque = cell(1,nL);
for j = 1:nL
    dL_dqd = diff(Lag, diff(th{j}(t), t));
    torque{j} = simplify(diff(dL_dqd, t) - diff(Lag, th{j}));
    fprintf('Torque %d derived\n', j);
end

%% Substitute g and rename kinematic variables
g_val = 9.81;
syms t1 t2 t3 t4 t5 t6 t1d t2d t3d t4d t5d t6d t1dd t2dd t3dd t4dd t5dd t6dd real
th_all      = [th1, th2, th3, th4, th5, th6];
th_dot_all  = [diff(th1,t), diff(th2,t), diff(th3,t), diff(th4,t), diff(th5,t), diff(th6,t)];
th_ddot_all = [diff(th1,t,t), diff(th2,t,t), diff(th3,t,t), diff(th4,t,t), diff(th5,t,t), diff(th6,t,t)];
for k = 1:nL
    tq = subs(torque{k}, g, g_val);
    tq = subs(tq, th_ddot_all, [t1dd t2dd t3dd t4dd t5dd t6dd]);   % 2nd derivatives first
    tq = subs(tq, th_dot_all,  [t1d t2d t3d t4d t5d t6d]);
    tq = subs(tq, th_all,      [t1 t2 t3 t4 t5 t6]);
    torque{k} = simplify(tq);
end

%% Save
for k = 1:nL
    fid = fopen(fullfile(out_folder, sprintf('Torque_%d.txt', k)), 'w');
    fprintf(fid, 'Torque_%d = %s;\n', k, char(torque{k}));
    fclose(fid);
    fprintf('Saved Torque_%d.txt\n', k);
end
