function V = exact_minimality(results_folder, prismatic, opts)
%EXACT_MINIMALITY Exact proof that the reduced regressor has full column rank, hence that p_b is minimal.
%
%   V = exact_minimality(results_folder)
%   V = exact_minimality(results_folder, prismatic)
%   V = exact_minimality(results_folder, prismatic, opts)
%
%   prismatic  optional n-by-1 logical, true for prismatic joints (the displacement of a prismatic joint
%              enters polynomially and is evaluated directly; only revolute angles use sin/cos).
%   opts       optional struct: seed (default 2024), n_states (default ceil(1.5*p_b/n)+4; doubled up to
%              twice if the rank is not yet full), control (default true, see below).
%
% Claim proved. Let A_I = A(:, independent_cols) be the reduced regressor (n-by-p_b), whose entries are
% polynomials in s_k = sin q_k, c_k = cos q_k, the velocities, the accelerations, sign(qd_k), the displacement
% of prismatic joints and, for the symbolic-dimension variants, the link dimensions. If the matrix obtained by
% stacking A_I at M exact states has rank p_b, the p_b columns are linearly independent as functions: no
% non-zero constant vector z satisfies A_I(x) z = 0 for all x, because it would also satisfy it at these
% states. Together with the exact dependency identity A = A_I*B (exact_verification) this gives
% rank A = p_b exactly, that is, p_b is the dimension of the identifiable subspace of the extracted regressor.
%
% Exactness. The states are exact rationals: each revolute angle is given by a Pythagorean pair
% (s, c) = (a/h, b/h) with a^2 + b^2 = h^2 and random signs, velocities, accelerations, displacements and
% dimensions are small rationals, and g_k = sign(qd_k) is the sign of the chosen velocity. Every entry of the
% stacked matrix is evaluated in exact rational arithmetic. Each row is multiplied by the least common
% multiple of its denominators (a non-zero factor, which does not change the rank) and reduced modulo the
% prime P = 2^26-5. Gauss elimination modulo P is exact in double precision because P < 2^26. A minor that is
% non-zero modulo P is non-zero over the rationals, hence rank_Q >= rank_P. A rank deficiency modulo P is not
% reported as a result: the test is repeated with more states, and if it persists the status is "not proved".
%
% Control. The complete regressor A (all p columns) is evaluated at the same states. Its rank must equal p_b
% as well, which checks the dependency certificate from the other side (rank A >= rank A_I = p_b and, by the
% certificate, rank A <= p_b).
% Writes exact_minimality.txt to results_folder.

    t0 = tic;
    S = load(fullfile(results_folder, 'base_regressor.mat'), 'A_matrix_common', 'independent_cols', 'dependent_cols');
    A = S.A_matrix_common;  ind = S.independent_cols(:).';  dep = S.dependent_cols(:).';
    [n, p] = size(A);  pb = numel(ind);
    if nargin < 2 || isempty(prismatic), prismatic = false(n, 1); end
    prismatic = logical(prismatic(:));
    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'seed'), opts.seed = 2024; end
    if ~isfield(opts, 'n_states'), opts.n_states = ceil(1.5*pb/n) + 4; end
    if ~isfield(opts, 'control'), opts.control = true; end
    P = 2^26 - 5;
    if ~isprime(P), error('exact_minimality:prime', 'P is not prime.'); end

    rev = find(~prismatic).';
    qn  = arrayfun(@(k) sym(sprintf('t%d', k), 'real'), 1:n);
    qdn = arrayfun(@(k) sym(sprintf('t%dd', k), 'real'), 1:n);
    s  = arrayfun(@(k) sym(sprintf('s%d', k), 'real'), rev);
    c  = arrayfun(@(k) sym(sprintf('c%d', k), 'real'), rev);
    sg = arrayfun(@(k) sym(sprintf('g%d', k), 'real'), 1:n);

    E = expand(A(:, [ind, dep]));                          % independent columns first
    E = subs(E, [sin(qn(rev)) cos(qn(rev)) sign(qdn)], [s c sg]);
    if any(has(E(:), [qn(rev), sign(qdn)]))
        error('exact_minimality:unsupported', 'Unsupported functions remain after expansion.');
    end
    vars = symvar(E);
    nv = numel(vars);
    [kind, joint] = classify(vars);

    pyth = [3 4 5; 5 12 13; 8 15 17; 7 24 25; 20 21 29; 12 35 37; 9 40 41; 28 45 53; 11 60 61; 16 63 65; 33 56 65; 48 55 73];
    rng(opts.seed, 'twister');
    M = opts.n_states;  attempt = 0;
    while true
        attempt = attempt + 1;
        Mrows = zeros(n*M, p);
        for m = 1:M
            qdv = randi(4, 1, n) .* (2*randi([0 1], 1, n) - 1);        % non-zero integers in [-4, 4]
            sv = sym(zeros(1, n));  cv = sym(zeros(1, n));
            for k = rev
                i = randi(size(pyth, 1));  a = pyth(i, 1);  b = pyth(i, 2);  h = pyth(i, 3);
                if randi(2) == 1, t = a; a = b; b = t; end
                sv(k) = sym(a*(2*randi([0 1]) - 1))/sym(h);
                cv(k) = sym(b*(2*randi([0 1]) - 1))/sym(h);
            end
            vals = sym(zeros(1, nv));
            for j = 1:nv
                k = joint(j);
                switch kind(j)
                    case 1, vals(j) = sv(k);                                   % s_k
                    case 2, vals(j) = cv(k);                                   % c_k
                    case 3, vals(j) = sym(sign(qdv(k)));                       % g_k = sign(qd_k)
                    case 4, vals(j) = sym(qdv(k));                             % qd_k
                    case 5, vals(j) = sym(randi([-4 4]));                      % qdd_k
                    case 6, vals(j) = sym(randi([-5 5]))/4;                    % displacement of a prismatic joint
                    otherwise, vals(j) = sym(randi([1 12]))/20;                % a link dimension
                end
            end
            Ev = subs(E, vars, vals);
            Mrows((m-1)*n + (1:n), :) = rows_mod_p(Ev, P);
        end
        rank_I = rank_mod_p(Mrows(:, 1:pb), P);
        if opts.control, rank_A = rank_mod_p(Mrows, P); else, rank_A = NaN; end
        if rank_I == pb || attempt >= 3, break; end
        M = 2*M;
    end
    V = struct('p', p, 'p_b', pb, 'n_states', M, 'n_rows', n*M, 'prime', P, 'rank_independent', rank_I, ...
        'rank_full', rank_A, 'proved', rank_I == pb, 'control_ok', ~opts.control || rank_A == pb, ...
        'attempts', attempt, 'seed', opts.seed, 'run_time_s', toc(t0));
    write_report(results_folder, V);
end

function [kind, joint] = classify(vars)
% kind: 1 s_k, 2 c_k, 3 g_k, 4 qd_k, 5 qdd_k, 6 displacement t_k, 7 any other symbol (a link dimension)
    nv = numel(vars);  kind = 7*ones(1, nv);  joint = zeros(1, nv);
    pat = {'^s(\d+)$', 1; '^c(\d+)$', 2; '^g(\d+)$', 3; '^t(\d+)d$', 4; '^t(\d+)dd$', 5; '^t(\d+)$', 6};
    for j = 1:nv
        nm = char(vars(j));
        for i = 1:size(pat, 1)
            tok = regexp(nm, pat{i, 1}, 'tokens', 'once');
            if ~isempty(tok), kind(j) = pat{i, 2};  joint(j) = str2double(tok{1});  break; end
        end
    end
end

function R = rows_mod_p(Ev, P)
% Rows of an exact rational matrix as integers modulo P (each row scaled by its common denominator).
    [N, D] = numden(Ev);
    R = zeros(size(Ev));
    for r = 1:size(Ev, 1)
        L = sym(1);
        for j = 1:size(Ev, 2), L = lcm(L, D(r, j)); end
        if mod(L, sym(P)) == 0, error('exact_minimality:prime', 'The prime divides a row scale; choose another prime.'); end
        row = N(r, :) .* (L ./ D(r, :));               % exact integers
        R(r, :) = double(mod(row, sym(P)));
    end
end

function r = rank_mod_p(M, P)
% Rank of an integer matrix modulo the prime P (< 2^26, so that products are exact in double precision).
    M = mod(M, P);  [m, nc] = size(M);  r = 0;
    for col = 1:nc
        if r == m, break; end
        piv = find(M(r+1:m, col) ~= 0, 1);
        if isempty(piv), continue; end
        piv = piv + r;
        M([r+1 piv], :) = M([piv r+1], :);
        r = r + 1;
        M(r, :) = mod(M(r, :) * inv_mod(M(r, col), P), P);
        rows = find(M(:, col) ~= 0);  rows(rows == r) = [];
        M(rows, :) = mod(M(rows, :) - M(rows, col) .* M(r, :), P);
    end
end

function x = inv_mod(a, P)
% a^(P-2) mod P by square and multiply; every product is below 2^52.
    x = 1;  e = P - 2;  b = mod(a, P);
    while e > 0
        if mod(e, 2) == 1, x = mod(x*b, P); end
        b = mod(b*b, P);  e = floor(e/2);
    end
end

function write_report(folder, V)
    fid = fopen(fullfile(folder, 'exact_minimality.txt'), 'w');
    fprintf(fid, 'Columns p / base columns p_b            : %d / %d\n', V.p, V.p_b);
    fprintf(fid, 'Exact states / stacked rows             : %d / %d (attempts %d, seed %d)\n', V.n_states, V.n_rows, V.attempts, V.seed);
    fprintf(fid, 'Prime for the exact elimination         : %d\n', V.prime);
    fprintf(fid, 'Rank of A_I at the states               : %d\n', V.rank_independent);
    fprintf(fid, 'Rank of the complete A at the states    : %d\n', V.rank_full);
    fprintf(fid, 'Full column rank of A_I proved exactly  : %d\n', V.proved);
    fprintf(fid, 'Control rank(A) = p_b                   : %d\n', V.control_ok);
    fprintf(fid, 'Run time                                : %.1f s\n', V.run_time_s);
    fclose(fid);
    type(fullfile(folder, 'exact_minimality.txt'));
end
