function [theta, A, info] = extract_regressor_columns_multi(tau, kin, dim_symbols, dim_values, known_flags, log_fn, use_parallel)
%EXTRACT_REGRESSOR_COLUMNS_MULTI Stage 1 for several dimension settings at once, optionally in parallel.
%
%   [theta, A, info] = extract_regressor_columns_multi(tau, kin, dim_symbols, dim_values, known_flags, log_fn, use_parallel)
%
%   known_flags  logical vector, one entry per setting (false: link dimensions
%                symbolic, true: substituted), for example [false true].
%   use_parallel true: the (setting, joint) pairs are distributed over the
%                workers of the current parallel pool with parfor. The pairs
%                are independent, and every one of them seeds its own random
%                states, so the result does not depend on the number of workers.
%
%   theta{s}{j}, A{s}{j}, info{s} follow extract_regressor_columns for setting s.
%   The method is that of extract_regressor_columns (derivative route).

    if nargin < 6 || isempty(log_fn), log_fn = @(varargin) []; end
    if nargin < 7, use_parallel = false; end
    n = numel(tau);  S = numel(known_flags);  T = n*S;
    res = cell(1, T);
    if use_parallel && T > 1
        parfor t = 1:T
            res{t} = joint_columns(tau(mod(t-1, n)+1), mod(t-1, n)+1, n, kin, dim_symbols, dim_values, ...
                known_flags(ceil(t/n)));
        end
    else
        for t = 1:T
            res{t} = joint_columns(tau(mod(t-1, n)+1), mod(t-1, n)+1, n, kin, dim_symbols, dim_values, ...
                known_flags(ceil(t/n)));
        end
    end
    theta = cell(1, S);  A = cell(1, S);  info = cell(1, S);
    for s = 1:S
        theta{s} = cell(1, n);  A{s} = cell(1, n);
        p_j = zeros(1, n);  chk = zeros(1, n);
        for j = 1:n
            r = res{(s-1)*n + j};
            theta{s}{j} = r.theta;  A{s}{j} = r.A;  p_j(j) = r.p;  chk(j) = r.chk;
            log_fn('  setting %d, joint %d: %d monomials, check %.1e\n', s, j, r.p, r.chk);
        end
        info{s} = struct('monomials_per_joint', p_j, 'check', chk);
    end
end

function r = joint_columns(Tj, j, n, kin, dim_symbols, dim_values, known)
% monomials and regressor columns of one torque equation (see extract_regressor_columns)
    if known && ~isempty(dim_symbols)
        Tj = subs(Tj, dim_symbols, sym(dim_values, 'r'));
    end
    phys = setdiff(symvar(Tj), kin, 'stable');
    mono = sym.empty(0, 1);
    for s = 1:2
        rng(100*j + s, 'twister');      % explicit generator: parallel-pool workers default to a different one
        state = [(-1+2*rand(1,n)) (-2+4*rand(1,n)) (-3+6*rand(1,n))];
        Te = expand(subs(Tj, kin, state));
        [~, mu] = coeffs(Te, phys);
        mono = [mono; mu(:)]; %#ok<AGROW>
    end
    [~, ui] = unique(arrayfun(@(e) string(char(e)), mono), 'stable');
    mono = mono(ui);
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
    theta = mono(keep);  Acol = cols(keep);
    rng(900 + j, 'twister');
    vals = [(-1+2*rand(1, 3*n)) 0.1+rand(1, numel(phys))];
    chk = abs(double(subs(Tj - Acol*theta, [kin phys], vals)));
    if chk > 1e-9
        error('extract_regressor_columns:check', 'Joint %d: columns do not reproduce the torque (%.3e).', j, chk);
    end
    r = struct('theta', {theta}, 'A', {Acol}, 'p', numel(theta), 'chk', chk);
end
