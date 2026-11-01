%% ========================================================================
% A2_regressor_columns_Panda.m
%
% Stage 1 of the regrouping for the Franka Emika Panda, carried out on the
% nested Euler-Lagrange expressions (panda_lagrange_torques.m) instead of the
% collected torque files.
%
% The collected Panda torques (torque_equations/Torque_k.txt, up to 4.6 MB
% per joint) are an order of magnitude larger than those of the UR5e and
% IRB 120; re-reading and expanding them dominates the run time. Since every
% torque is a polynomial in the physical symbols, the regressor column of a
% monomial mu = prod x^alpha is obtained directly as
%       c(q,qd,qdd) = (1/alpha!) d^alpha tau / dx^alpha  at x = 0,
% which is identical to the column of the placeholder Jacobian used for the
% other robots. Monomials are found exactly as in Stage 1 (coefficients of
% tau at numerical kinematic states).
%
% Two cases are written: link dimensions symbolic, and known (numeric).
% OUTPUT: torque_equations/regressor_columns.mat
%         torque_equations/regressor_columns_known_geometry.mat
%         (theta{j}: monomials of joint j, A{j}: their columns)
% Run time: see A2_log.txt.
%% ========================================================================

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here);
out_folder = fullfile(here, 'torque_equations');
log_file = fullfile(here, 'A2_log.txt');
if isfile(log_file), delete(log_file); end
t_start = tic;

nL = 7;
q   = arrayfun(@(k) sym(sprintf('t%d',   k), 'real'), 1:nL);
qd  = arrayfun(@(k) sym(sprintf('t%dd',  k), 'real'), 1:nL);
qdd = arrayfun(@(k) sym(sprintf('t%ddd', k), 'real'), 1:nL);
syms L1 L2 L3 L4 L5 L6 real
geometry_symbols = [L1 L2 L3 L4 L5 L6];
geometry_values  = sym([333 316 82.5 82.5 384 88])/1000;
g_val = sym(981)/100;

%% Euler-Lagrange equations in the parallel-axis form (as in script A)
tau_nested = panda_lagrange_torques(q, qd, qdd, geometry_symbols, g_val);
logp(log_file, 'Euler-Lagrange equations (%.0f s)\n', toc(t_start));

%% Monomials and regressor columns, symbolic and known geometry
kin = [q qd qdd];
cases = {'regressor_columns', false; 'regressor_columns_known_geometry', true};
for c = 1:size(cases, 1)
    theta = cell(1, nL);  A = cell(1, nL);
    for j = 1:nL
        Tj = tau_nested(j);
        if cases{c, 2}, Tj = subs(Tj, geometry_symbols, geometry_values); end
        phys = setdiff(symvar(Tj), kin, 'stable');
        % monomials: coefficients of tau at two numerical kinematic states
        mono = sym.empty(0,1);
        for s = 1:2
            rng(100*j + s);
            Te = expand(subs(Tj, kin, [(-1+2*rand(1,nL)) (-2+4*rand(1,nL)) (-3+6*rand(1,nL))]));
            [~, mu] = coeffs(Te, phys);
            mono = [mono; mu(:)]; %#ok<AGROW>
        end
        [~, ui] = unique(arrayfun(@(e) string(char(e)), mono), 'stable');
        mono = mono(ui);
        % column of each monomial: (1/alpha!) d^alpha tau/dx^alpha at x = 0
        zero_phys = zeros(1, numel(phys));
        cols = sym(zeros(1, numel(mono)));
        for k = 1:numel(mono)
            D = Tj;  fact = 1;
            for v = 1:numel(phys)
                e = polynomialDegree(mono(k), phys(v));
                if e > 0, D = diff(D, phys(v), e);  fact = fact*factorial(e); end
            end
            cols(k) = expand(subs(D, phys, zero_phys)/fact);
        end
        keep = arrayfun(@(e) ~isequal(e, sym(0)), cols);
        theta{j} = mono(keep);  A{j} = cols(keep);
        % the columns must reproduce the nested Euler-Lagrange torque
        rng(900 + j);
        vals = [(-1+2*rand(1,3*nL)) 0.1+rand(1,numel(phys))];
        e_chk = double(subs(Tj - A{j}*theta{j}, [kin phys], vals));
        if abs(e_chk) > 1e-9, error('Joint %d: columns do not reproduce the torque (%.3e).', j, e_chk); end
        logp(log_file, '%s, joint %d: %d monomials, check %.1e (%.0f s)\n', ...
            cases{c, 1}, j, numel(theta{j}), abs(e_chk), toc(t_start));
    end
    save(fullfile(out_folder, [cases{c, 1} '.mat']), 'theta', 'A');
end
logp(log_file, 'Total run time %.0f s\n', toc(t_start));

function logp(file, varargin)
    % progress is appended to the log file line by line
    fid = fopen(file, 'a');  fprintf(fid, varargin{:});  fclose(fid);
    fprintf(varargin{:});
end
