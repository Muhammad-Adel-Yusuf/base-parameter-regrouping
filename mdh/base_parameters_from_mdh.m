function out = base_parameters_from_mdh(mdh, opts)
%BASE_PARAMETERS_FROM_MDH Certified base parameters of a serial manipulator from its MDH table.
%
%   out = base_parameters_from_mdh(mdh)
%   out = base_parameters_from_mdh(mdh, opts)
%
%   mdh   struct from mdh_table.m (modified DH table, revolute and prismatic
%         joints; the inertial parameters stay symbolic)
%   opts  optional struct:
%     results_root          folder for the results (default: ../results)
%     variants              cell array, any of
%                             'symbolic'         link dimensions symbolic
%                             'known'            link dimensions substituted
%                             'actuation'        symbolic, with Ia, Fv, Fc
%                             'actuation_known'  known, with Ia, Fv, Fc
%                           (default: all four; the known variants need dims)
%     exact_reconstruction  true: also prove tau_sym = A*Theta exactly by
%                           polynomial reduction (default: n <= 6; it expands
%                           the full torque and is slow for large models)
%     exact_minimality      true: prove exactly that the p_b reduced columns are
%                           linearly independent (exact_minimality.m; default true)
%     parallel              true: use a parallel pool (default true when the
%                           Parallel Computing Toolbox is available): Stage 1 over
%                           (setting, joint) and the variants side by side
%     parallel_variants     true: run the variants in parallel (default: n <= 6, or
%                           at least 48 GB of free memory; the four variants of the
%                           7-joint Panda peaked at 33.5 GB in total)
%     workers               pool size if a pool has to be started (default 10)
%     cache                 true: reuse the Stage 1 columns of an earlier run of
%                           the same table (default true); the cache is the
%                           regressor_columns*.mat file with a key of the table
%
%   Steps (all automatic):
%     1  symbolic Euler-Lagrange torques from the table      mdh_lagrange_torques
%     2  Stage 1-2: monomials and regressor columns          extract_regressor_columns_multi
%     3  Stages 3-5: rank, pivoted QR, regrouping, numerical
%        verification against an independent RNEA            regroup_base_parameters
%     4  exact dependency certificate (polynomial reduction) exact_verification
%     5  exact reconstruction certificate (optional)
%     6  exact minimality proof: full column rank of A_b at exact rational states
%     7  base count with numerical dimensions, once more from the RNEA alone
%
%   Returns out.table (one row per variant), out.folder and out.run_time_s.
%   Writes <results_root>/<name>/<variant>/ with theta_base_list.txt,
%   regrouping_summary.txt, exact_verification.txt and base_regressor.mat, and
%   <results_root>/<name>/summary.txt.

    here = fileparts(mfilename('fullpath'));
    core = fullfile(fileparts(here), 'core');
    addpath(here);  addpath(core);
    if nargin < 2, opts = struct(); end
    def = struct('results_root', fullfile(fileparts(here), 'results'), ...
        'variants', {{'symbolic', 'known', 'actuation', 'actuation_known'}}, ...
        'exact_reconstruction', mdh.n <= 6, 'exact_minimality', true, 'parallel', true, 'parallel_variants', mdh.n <= 6 || enough_memory(48e9), ...
        'workers', 10, 'cache', true);
    fn = fieldnames(def);
    for k = 1:numel(fn), if ~isfield(opts, fn{k}), opts.(fn{k}) = def.(fn{k}); end, end

    n = mdh.n;
    folder = fullfile(opts.results_root, mdh.name);
    if ~isfolder(folder), mkdir(folder); end
    log_file = fullfile(folder, 'run_log.txt');
    if isfile(log_file), delete(log_file); end
    logp = @(varargin) logline(log_file, varargin{:});
    t_all = tic;

    use_par = opts.parallel && start_pool(opts.workers, here, core, logp);

    has_dims = ~isempty(mdh.dim_names);
    variants = opts.variants;
    if ~has_dims
        variants = variants(~ismember(variants, {'known', 'actuation_known'}));
        logp('No link dimensions declared: known-geometry variants skipped.\n');
    end
    folder_of = struct('symbolic', 'symbolic', 'known', 'known_geometry', ...
        'actuation', 'actuation', 'actuation_known', 'actuation_known_geometry');
    known_of  = struct('symbolic', false, 'known', true, 'actuation', false, 'actuation_known', true);
    act_of    = struct('symbolic', false, 'known', false, 'actuation', true, 'actuation_known', true);

    %% 1  symbolic torques from the table
    logp('[%s] %d joints (%d prismatic)\n', mdh.name, n, sum(mdh.prismatic));
    t1 = tic;
    [tau, q, qd, qdd] = mdh_lagrange_torques(mdh);
    kin = [q qd qdd];
    t_torque = toc(t1);
    logp('Euler-Lagrange torques derived (%.1f s)\n', t_torque);

    %% 2  Stage 1-2 extraction (once per dimension setting; cached per table)
    need_known = any(cellfun(@(v) known_of.(v), variants));
    need_sym   = any(cellfun(@(v) ~known_of.(v), variants));
    key = table_key(mdh);
    colfile = struct('symbolic', fullfile(folder, 'regressor_columns.mat'), ...
                     'known',    fullfile(folder, 'regressor_columns_known_geometry.mat'));
    flags = false(1, 0);  files = {};
    if need_sym,   flags(end+1) = false;  files{end+1} = colfile.symbolic; end %#ok<AGROW>
    if need_known, flags(end+1) = true;   files{end+1} = colfile.known;    end %#ok<AGROW>
    todo = true(1, numel(files));
    if opts.cache
        for k = 1:numel(files)
            todo(k) = ~cache_valid(files{k}, key);
            if ~todo(k), logp('Stage 1-2 (%s dimensions): reused from the cache\n', ternary(flags(k), 'known', 'symbolic')); end
        end
    end
    t2 = tic;
    if any(todo)
        logp('Stage 1-2 for %s dimensions (%s)\n', strjoin(arrayfun(@(f) ternary(f, 'known', 'symbolic'), flags(todo), 'UniformOutput', false), ' and '), ...
            ternary(use_par, 'parallel', 'sequential'));
        [th, AA] = extract_regressor_columns_multi(tau, kin, mdh.dim_symbols, mdh.dim_values, flags(todo), logp, use_par);
        idx = find(todo);
        for k = 1:numel(idx)
            theta = th{k};  A = AA{k};                                       %#ok<NASGU>
            save(files{idx(k)}, 'theta', 'A', 'key');
        end
    end
    t_extract = toc(t2);

    %% 3-5  regrouping, dependency certificate, reconstruction certificate (+ independent count)
    if has_dims
        geom_names = mdh.dim_names;  geom_values = mdh.dim_values;
    else
        geom_names = {'Ldummy'};  geom_values = 0;
    end
    r_rnea = struct('rigid', NaN, 'act', NaN);          % independent count, numerical dimensions
    if has_dims
        r_rnea.rigid = mdh_numeric_rank(mdh, false);
        r_rnea.act   = mdh_numeric_rank(mdh, true);
    end
    nv = numel(variants);
    rows_c = cell(1, nv);  msgs = cell(1, nv);
    ctx = struct('mdh', mdh, 'tau', {tau}, 'q', {q}, 'qd', {qd}, 'qdd', {qdd}, 'folder', folder, ...
        'folder_of', folder_of, 'known_of', known_of, 'act_of', act_of, 'colfile', colfile, ...
        'geom_names', {geom_names}, 'geom_values', geom_values, 'r_rnea', r_rnea, ...
        'exact_reconstruction', opts.exact_reconstruction, 'exact_minimality', opts.exact_minimality);
    if use_par && opts.parallel_variants && nv > 1
        parfor iv = 1:nv
            [row, msg] = run_variant(variants{iv}, ctx);
            rows_c{iv} = row;  msgs{iv} = msg;
        end
    else
        for iv = 1:nv
            [row, msg] = run_variant(variants{iv}, ctx);
            rows_c{iv} = row;  msgs{iv} = msg;
        end
    end
    rows = [rows_c{:}];
    for iv = 1:nv, logp('\n== variant %s ==\n%s', variants{iv}, msgs{iv}); end

    %% summary
    out = struct('name', mdh.name, 'folder', folder, 'table', rows, ...
        'torque_derivation_s', t_torque, 'extraction_s', t_extract, 'run_time_s', toc(t_all), ...
        'parallel', use_par);
    save(fullfile(folder, 'summary.mat'), 'out');
    fid = fopen(fullfile(folder, 'summary.txt'), 'w');
    fprintf(fid, '%s (%d joints, %d prismatic)\n', mdh.name, n, sum(mdh.prismatic));
    fprintf(fid, '%-18s %6s %6s %12s %16s %12s %12s %8s %9s %10s\n', 'variant', 'p', 'p_b', 'dependency', 'reconstruction', ...
        'minimality', 'e_max (Nm)', 'r_geom', 'r_RNEA', 'time (s)');
    for k = 1:numel(rows)
        fprintf(fid, '%-18s %6d %6d %12s %16s %12s %12.2e %8d %9d %10.0f\n', rows(k).variant, rows(k).p, rows(k).p_b, ...
            ternary(rows(k).dependency_certified, 'certified', 'NO'), rows(k).reconstruction, rows(k).minimality, ...
            rows(k).e_max, rows(k).r_geom, rows(k).r_geom_rnea, rows(k).time_s);
    end
    fprintf(fid, ['r_geom = base count with the link dimensions substituted (from the regrouping); r_RNEA = the same count from the\n' ...
        'Newton-Euler reference alone (mdh_numeric_rank), which does not use the symbolic torques.\n']);
    fprintf(fid, 'torque derivation %.1f s; extraction %.1f s (%s); total %.0f s (%s)\n', ...
        t_torque, t_extract, ternary(any(todo), 'computed', 'from cache'), out.run_time_s, ...
        ternary(use_par, 'parallel', 'sequential'));
    fclose(fid);
    logp('\nTotal run time %.0f s\n', out.run_time_s);
end

function [row, msg] = run_variant(v, c)
% One variant: regrouping, dependency certificate, reconstruction certificate. Returns the table row
% and its log text (no file logging here, so it can run inside parfor).
    tv = tic;
    mdh = c.mdh;
    ref = mdh_rnea_reference(mdh);
    rf = fullfile(c.folder, c.folder_of.(v));
    opt = struct('robot_name', mdh.name, 'nDOF', mdh.n, 'torque_folder', '', 'results_folder', rf, ...
        'geometry_names', {c.geom_names}, 'geometry_values', c.geom_values, ...
        'known_geometry', c.known_of.(v), 'actuation', c.act_of.(v), ...
        'torque_function', ref);
    if c.known_of.(v), opt.column_file = c.colfile.known; else, opt.column_file = c.colfile.symbolic; end
    s = regroup_base_parameters(opt);
    V = exact_verification(rf, mdh.prismatic);
    recon = 'not run';  t_recon = 0;
    if c.exact_reconstruction
        tr = tic;
        recon = ternary(exact_reconstruction(c.tau, c.q, c.qd, c.qdd, mdh, c.known_of.(v), c.act_of.(v), rf), ...
            'certified', 'NOT certified');
        t_recon = toc(tr);
    end
    mini = 'not run';  t_min = 0;
    if c.exact_minimality
        Mi = exact_minimality(rf, mdh.prismatic);
        t_min = Mi.run_time_s;
        if ~Mi.proved || ~Mi.control_ok
            mini = 'NOT proved';
        elseif V.certified || s.p_b == s.p
            mini = 'proved';                       % lower bound (independence) and upper bound (dependency identity)
        else
            mini = 'lower bound';                  % independence proved, dependency identity not certified
        end
    end
    r_rnea = ternary(c.act_of.(v), c.r_rnea.act, c.r_rnea.rigid);
    row = struct('variant', v, 'p', s.p, 'p_b', s.p_b, ...
        'dependency_certified', V.certified, 'reconstruction', recon, 'minimality', mini, ...
        'e_max', s.e_max, 'e_full', s.e_full, 'r_geom', s.r_numeric_geometry, ...
        'r_geom_rnea', r_rnea, 'stage_times_s', s.stage_times_s, 'dependency_s', V.run_time_s, ...
        'reconstruction_s', t_recon, 'minimality_s', t_min, 'time_s', toc(tv));
    msg = sprintf('variant %s: p = %d, p_b = %d, dependency certified = %d, reconstruction %s, minimality %s, numeric-geometry count %d (RNEA alone: %d) (%.0f s)\n', ...
        v, s.p, s.p_b, V.certified, recon, mini, s.r_numeric_geometry, r_rnea, row.time_s);
end

function ok = exact_reconstruction(tau, q, qd, qdd, mdh, known, act, rfolder)
% Exact proof of tau_sym - A*Theta = 0 by reduction modulo s_k^2 + c_k^2 - 1.
    n = mdh.n;
    rev = find(~mdh.prismatic(:)).';              % prismatic displacements stay polynomial variables
    S = load(fullfile(rfolder, 'base_regressor.mat'), 'A_matrix_common', 'theta_common');
    A = S.A_matrix_common;  th = S.theta_common(:);
    s  = arrayfun(@(k) sym(sprintf('s%d', k), 'real'), rev);
    c  = arrayfun(@(k) sym(sprintf('c%d', k), 'real'), rev);
    sg = arrayfun(@(k) sym(sprintf('g%d', k), 'real'), 1:n);
    ideal = s.^2 + c.^2 - 1;
    ok = true;
    for j = 1:n
        Tj = tau(j);
        if known && ~isempty(mdh.dim_symbols), Tj = subs(Tj, mdh.dim_symbols, sym(mdh.dim_values, 'r')); end
        if act
            Tj = Tj + sym(sprintf('Ia%d', j), 'real')*qdd(j) + sym(sprintf('Fv%d', j), 'real')*qd(j) ...
                    + sym(sprintf('Fc%d', j), 'real')*sign(qd(j));
        end
        Rj = expand(Tj - A(j, :)*th);
        Rj = subs(Rj, [sin(q(rev)) cos(q(rev)) sign(qd)], [s c sg]);
        if any(has(Rj, [q(rev), sign(qd)])), ok = false; return; end
        if isequal(Rj, sym(0)), continue; end
        vars = [s c setdiff(symvar(Rj), [s c])];
        r = polynomialReduce(Rj, ideal, vars, 'MonomialOrder', 'lexicographic');
        if ~isequal(expand(r), sym(0)), ok = false; return; end
    end
    fid = fopen(fullfile(rfolder, 'reconstruction_check.txt'), 'w');
    fprintf(fid, 'Joint residuals tau_sym - A*Theta tested : %d\nAll remainders zero                      : 1\n', n);
    fclose(fid);
end

function ok = start_pool(workers, here, core, logp)
% Make sure a process pool with the project folders on the workers' path is running.
    ok = false;
    try
        if ~license('test', 'Distrib_Computing_Toolbox'), logp('Parallel Computing Toolbox not available: sequential run.\n'); return; end
        pool = gcp('nocreate');
        if isempty(pool)
            workers = min(workers, feature('numcores'));
            pool = parpool('Processes', workers);
        end
        pctRunOnAll(sprintf('addpath(''%s'', ''%s'')', here, core));
        logp('Parallel pool with %d workers.\n', pool.NumWorkers);
        ok = pool.NumWorkers > 1;
    catch err
        logp('Parallel pool could not be started (%s): sequential run.\n', err.message);
    end
end

function ok = enough_memory(bytes)
% true if at least `bytes` of physical memory are available (Windows); false where this cannot be told
    try
        [~, sys] = memory;
        ok = sys.PhysicalMemory.Available >= bytes;
    catch
        ok = false;
    end
end

function key = table_key(mdh)
% identifies the table (and the extraction method) for the Stage 1 cache
    items = [mdh.alpha(:); mdh.a(:); mdh.d(:); mdh.theta0(:)].';
    s = [mdh.name '|' strjoin(arrayfun(@char, items, 'UniformOutput', false), ',') '|' mat2str(double(mdh.prismatic(:).')) ...
        '|' strjoin(mdh.dim_names, ',') '|' mat2str(mdh.dim_values, 17) '|extract-v1'];
    try
        md = java.security.MessageDigest.getInstance('MD5');
        key = lower(reshape(dec2hex(typecast(md.digest(uint8(s)), 'uint8'), 2).', 1, []));
    catch
        key = sprintf('%d-%.0f', numel(s), mod(sum(double(s).*(1:numel(s))), 2^50));
    end
end

function ok = cache_valid(file, key)
    ok = false;
    if ~isfile(file), return; end
    w = whos('-file', file);
    if ~any(strcmp({w.name}, 'key')), return; end
    S = load(file, 'key');
    ok = strcmp(S.key, key);
end

function logline(file, varargin)
    msg = sprintf(varargin{:});
    fid = fopen(file, 'a');  fprintf(fid, '%s', msg);  fclose(fid);
    fprintf('%s', msg);
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end
