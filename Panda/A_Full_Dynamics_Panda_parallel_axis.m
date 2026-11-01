%% ========================================================================
% A_Full_Dynamics_Panda_parallel_axis.m
%
% Symbolic Euler-Lagrange inverse dynamics of the Franka Emika Panda (7 DOF,
% modified DH, Craig convention as published by Franka Emika) using
% ORIGIN-REFERENCED inertia tensors (parallel-axis formulation), as for the
% UR5e and the IRB 120:
%
%   K_i = 1/2 m_i v_i'v_i + m_i v_i'(w_i x c_i) + 1/2 w_i' I_oi w_i
%   V_i = -m_i g0'(R_i c_i + o_i)
%   tau_j = d/dt(dL/dqd_j) - dL/dq_j
%
% where I_oi = I_ci + m_i (c_i'c_i E3 - c_i c_i') is the inertia tensor of
% link i about the ORIGIN of frame {i} (symbols Ixx_i..Izz_i). The torques
% are then linear in the 10 standard parameters per link. The derivation is
% carried out by panda_lagrange_torques.m without symbolic simplification
% (seconds); the symbolic simplification used for the 6-DOF robots is not
% feasible for the Panda, whose collected torques are an order of magnitude
% larger.
%
% Each torque is then written in collected form, sum_k c_k(q,qd,qdd)*mu_k,
% where mu_k are products of the physical symbols (m_i, Lc_i, Ixx_i.., L1..L6)
% and the trigonometric coefficients c_k are expanded one by one into
% polynomials in the sines and cosines of the joint angles. The coefficient
% of mu_k = prod x^alpha is obtained without expanding the nested expression
% as (1/alpha!) d^alpha tau / dx^alpha at x = 0, and the collected form is
% checked against the nested one at random values.
%
% g is substituted numerically (known constant). Link dimensions L1..L6 are
% kept symbolic: d1 = L1, d3 = L2, a3 = L3, a4 = -L4, d5 = L5, a6 = L6.
%
% OUTPUT: torque_equations/Torque_1.txt ... Torque_7.txt
% Run time: see A_log.txt (coefficient expansion dominates).
%% ========================================================================

clear; clc;

here = fileparts(mfilename('fullpath'));
addpath(here);
out_folder = fullfile(here, 'torque_equations');
if ~isfolder(out_folder), mkdir(out_folder); end
log_file = fullfile(here, 'A_log.txt');
if isfile(log_file), delete(log_file); end
t_start = tic;

nL = 7;
q   = arrayfun(@(k) sym(sprintf('t%d',   k), 'real'), 1:nL);
qd  = arrayfun(@(k) sym(sprintf('t%dd',  k), 'real'), 1:nL);
qdd = arrayfun(@(k) sym(sprintf('t%ddd', k), 'real'), 1:nL);
syms L1 L2 L3 L4 L5 L6 real
g_val = sym(981)/100;

%% Euler-Lagrange equations (energies in the parallel-axis form)
tau_nested = panda_lagrange_torques(q, qd, qdd, [L1 L2 L3 L4 L5 L6], g_val);
logp(log_file, 'Euler-Lagrange equations (%.0f s)\n', toc(t_start));

%% Collected form: tau_j = sum_k c_k * mu_k
kin = [q qd qdd];
torque = cell(1,nL);
for j = 1:nL
    Tj = tau_nested(j);
    phys = setdiff(symvar(Tj), kin, 'stable');
    % monomials from two random kinematic states (a coefficient vanishing at
    % both is a measure-zero event; the check below would detect it)
    mono = sym.empty(0,1);
    for s = 1:2
        rng(100*j + s);
        Te = expand(subs(Tj, kin, [(-1+2*rand(1,nL)) (-2+4*rand(1,nL)) (-3+6*rand(1,nL))]));
        [~, mu] = coeffs(Te, phys);
        mono = [mono; mu(:)]; %#ok<AGROW>
    end
    [~, ui] = unique(arrayfun(@(e) string(char(e)), mono), 'stable');
    mono = mono(ui);
    zero_phys = zeros(1, numel(phys));
    Tc = sym(0);
    for k = 1:numel(mono)
        D = Tj;  fact = 1;
        for v = 1:numel(phys)
            e = polynomialDegree(mono(k), phys(v));
            if e > 0, D = diff(D, phys(v), e);  fact = fact*factorial(e); end
        end
        ck = expand(subs(D, phys, zero_phys)/fact);
        Tc = Tc + ck*mono(k);
    end
    % check the collected form against the nested Euler-Lagrange expression
    rng(900 + j);
    vals = [(-1+2*rand(1,3*nL)) 0.1+rand(1,numel(phys))];
    e_chk = double(subs(Tj - Tc, [kin phys], vals));
    if abs(e_chk) > 1e-9, error('Torque %d: collected form differs (%.3e).', j, e_chk); end
    torque{j} = Tc;
    logp(log_file, 'Torque %d: %d monomials, check %.1e (%.0f s)\n', j, numel(mono), abs(e_chk), toc(t_start));
end

%% Save
for k = 1:nL
    fid = fopen(fullfile(out_folder, sprintf('Torque_%d.txt', k)), 'w');
    fprintf(fid, 'Torque_%d = %s;\n', k, char(torque{k}));
    fclose(fid);
    logp(log_file, 'Saved Torque_%d.txt\n', k);
end
logp(log_file, 'Total run time %.0f s\n', toc(t_start));

function logp(file, varargin)
    % progress is appended to the log file line by line
    fid = fopen(file, 'a');  fprintf(fid, varargin{:});  fclose(fid);
    fprintf(varargin{:});
end
