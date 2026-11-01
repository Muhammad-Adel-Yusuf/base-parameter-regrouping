function summary = regroup_base_parameters(opt)
%REGROUP_BASE_PARAMETERS Automatic base-parameter regrouping (Algorithm 1).
%
%   summary = regroup_base_parameters(opt)
%
%   Stage 1  symbolic monomial extraction from each torque equation
%   Stage 2  common multi-joint regressor A(x) (6 x p)
%   Stage 3  structural rank from M random states
%   Stage 4  column-pivoted QR -> A_b(x), Theta_b = Theta_I + Gamma*Theta_D
%   Stage 5  verification at independent states (random physical values)
%   Extra    base count with the link dimensions substituted numerically
%
%   opt fields: robot_name, torque_folder, results_folder, geometry_names,
%   geometry_values, known_geometry (true: substitute the link dimensions
%   before extraction, i.e. regroup with known geometry), and optionally
%   nDOF (6), discovery_seed_offset (10), rank_seed (777),
%   verification_seed (4242), n_verification_states (200), rank_factor (1e-8),
%   actuation (false: rigid-body torques only; true: add the rotor inertia,
%   viscous and Coulomb friction of every joint, Ia_j*qdd_j + Fv_j*qd_j +
%   Fc_j*sign(qd_j), as symbolic parameters).

    def = struct('nDOF', 6, 'discovery_seed_offset', 10, 'rank_seed', 777, ...
        'verification_seed', 4242, 'n_verification_states', 200, ...
        'rank_factor', 1e-8, 'known_geometry', false, 'actuation', false);
    fn = fieldnames(def);
    for k = 1:numel(fn)
        if ~isfield(opt, fn{k}), opt.(fn{k}) = def.(fn{k}); end
    end
    robot_name = opt.robot_name;  torque_folder = opt.torque_folder;
    results_folder = opt.results_folder;  nDOF = opt.nDOF;
    geometry_names = opt.geometry_names;  geometry_values = opt.geometry_values;
    discovery_seed_offset = opt.discovery_seed_offset;  rank_seed = opt.rank_seed;
    verification_seed = opt.verification_seed;
    n_verification_states = opt.n_verification_states;  rank_factor = opt.rank_factor;

    if ~isfolder(results_folder), mkdir(results_folder); end
    t_start = tic;

    syms t1 t2 t3 t4 t5 t6 t1d t2d t3d t4d t5d t6d t1dd t2dd t3dd t4dd t5dd t6dd real
    q   = [t1 t2 t3 t4 t5 t6];
    qd  = [t1d t2d t3d t4d t5d t6d];
    qdd = [t1dd t2dd t3dd t4dd t5dd t6dd];
    kinematic_variables = [q qd qdd];

    %% ========================================================================
    % STAGE 1-2: MONOMIAL EXTRACTION AND JOINT REGRESSORS
    % ========================================================================
    torques = load_torque_equations(torque_folder, nDOF);
    if opt.known_geometry
        % regroup with known link dimensions (classical setting)
        for jn = 1:nDOF
            torques{jn} = subs(torques{jn}, sym(geometry_names), sym(geometry_values));
        end
    end
    actuation_symbols = sym([]);
    if opt.actuation
        % joint-side actuation terms: rotor inertia, viscous and Coulomb friction
        for jn = 1:nDOF
            Ia = sym(sprintf('Ia%d', jn), 'real');
            Fv = sym(sprintf('Fv%d', jn), 'real');
            Fc = sym(sprintf('Fc%d', jn), 'real');
            torques{jn} = torques{jn} + Ia*qdd(jn) + Fv*qd(jn) + Fc*sign(qd(jn));
            actuation_symbols = [actuation_symbols, Ia, Fv, Fc]; %#ok<AGROW>
        end
    end

    results = struct('Torque', {}, 'theta', {}, 'A', {});
    p_i = zeros(1, nDOF);
    for jn = 1:nDOF
        Ti = torques{jn};
        physical_symbols = setdiff(symvar(Ti), kinematic_variables, 'stable');

        % discovery state: kinematics numeric, physical symbols kept symbolic
        rng(discovery_seed_offset + jn);
        kv = [(-1+2*rand(1,nDOF)) (-2+4*rand(1,nDOF)) (-3+6*rand(1,nDOF))];
        Te = expand(subs(Ti, kinematic_variables, kv));
        if isequal(Te, sym(0))
            term_syms = sym.empty(0,1);
        elseif isequal(char(feval(symengine, 'op', Te, 0)), '_plus')
            term_syms = children(Te);
            term_syms = reshape([term_syms{:}], [], 1);
        else
            term_syms = Te;
        end

        % separate the numerical coefficient and the physical monomial of each
        % term; coeffs() returns the monomial as a pure product of physical
        % symbols, so that numerical factors left unevaluated by the
        % substitution (e.g. products of sines of the discovery state) cannot
        % create spurious duplicates of the same monomial
        parameter_list = sym.empty(0,1);
        for k = 1:numel(term_syms)
            [c_k, mono] = coeffs(term_syms(k), physical_symbols);
            if isempty(mono) || isequal(mono(1), sym(1)) || isAlways(c_k(1) == 0, 'Unknown', 'false')
                continue
            end
            parameter_list(end+1,1) = mono(1); %#ok<AGROW>
        end
        [~, ui] = unique(arrayfun(@(e) string(char(e)), parameter_list), 'stable');
        theta_i = parameter_list(ui);
        np = numel(theta_i);

        % placeholder substitution (longest expression first) and Jacobian
        P = sym('P', [np 1], 'real');
        [~, ord] = sort(arrayfun(@(e) length(char(e)), theta_i), 'descend');
        Torque_P = subs(expand(Ti), theta_i(ord), P(ord));
        if ~isempty(intersect(symvar(Torque_P), physical_symbols))
            for k = 1:np
                if ~isempty(intersect(symvar(Torque_P), symvar(theta_i(ord(k)))))
                    Torque_P = subs(Torque_P, theta_i(ord(k)), P(ord(k)));
                end
            end
            leftover = intersect(symvar(Torque_P), physical_symbols);
            if ~isempty(leftover)
                error('Joint %d: could not substitute %s.', jn, ...
                    strjoin(arrayfun(@char, leftover, 'UniformOutput', false), ', '));
            end
        end
        results(jn).Torque = Ti;
        results(jn).theta  = theta_i;
        results(jn).A      = jacobian(Torque_P, P);
        p_i(jn) = np;
        fprintf('Joint %d: %d monomials (%.0f s)\n', jn, np, toc(t_start));
    end

    %% Common parameter vector and aligned regressor
    all_theta = vertcat(results.theta);
    all_strings = arrayfun(@(e) string(char(e)), all_theta);
    [~, ui] = unique(all_strings, 'stable');
    theta_common = all_theta(ui);
    theta_common_strings = all_strings(ui);
    p = numel(theta_common);

    A_matrix_common = sym(zeros(nDOF, p));
    for jn = 1:nDOF
        s = arrayfun(@(e) string(char(e)), results(jn).theta);
        for k = 1:numel(s)
            A_matrix_common(jn, find(theta_common_strings == s(k), 1)) = results(jn).A(k);
        end
    end
    A_common_fn = matlabFunction(A_matrix_common, 'Vars', {q, qd, qdd});
    t_stage12 = toc(t_start);

    %% ========================================================================
    % STAGE 3: STRUCTURAL RANK
    % ========================================================================
    M = max(50, ceil(10*p/nDOF));
    Phi_str = zeros(nDOF*M, p);
    rng(rank_seed);
    for j = 1:M
        qj = -pi + 2*pi*rand(1,nDOF); qdj = -2 + 4*rand(1,nDOF); qddj = -3 + 6*rand(1,nDOF);
        Phi_str(nDOF*(j-1)+(1:nDOF), :) = A_common_fn(qj, qdj, qddj);
    end
    sv_str = svd(Phi_str);
    rank_tol = rank_factor * max(size(Phi_str)) * max(sv_str);
    r = sum(sv_str > rank_tol);
    fprintf('\np = %d, M = %d, structural rank r = %d\n', p, M, r);
    t_stage3 = toc(t_start) - t_stage12;

    %% ========================================================================
    % STAGE 4: COLUMN-PIVOTED QR AND SYMBOLIC REGROUPING
    % ========================================================================
    if r == p
        theta_base = theta_common;  A_matrix_base = A_matrix_common;
        independent_cols = 1:p;  dependent_cols = [];  Gamma = zeros(p,0);
    else
        [~, R, pivot] = qr(Phi_str, 0);
        independent_cols = pivot(1:r);
        dependent_cols   = pivot(r+1:end);
        Gamma = R(1:r,1:r) \ R(1:r,r+1:end);
        Gamma(abs(Gamma) < 1e4*eps(class(Gamma))*max(abs(Gamma(:)))) = 0;
        theta_base    = theta_common(independent_cols) + Gamma*theta_common(dependent_cols);
        A_matrix_base = A_matrix_common(:, independent_cols);
    end
    p_b = numel(theta_base);
    A_base_fn = matlabFunction(A_matrix_base, 'Vars', {q, qd, qdd});
    t_stage4 = toc(t_start) - t_stage12 - t_stage3;

    %% ========================================================================
    % STAGE 5: VERIFICATION AT INDEPENDENT STATES
    % ========================================================================
    physical_all = symvar(theta_common);
    rng(verification_seed);
    pv = 0.1 + rand(1, numel(physical_all));               % random physical values
    theta_common_num = double(subs(theta_common, physical_all, pv));
    theta_base_num   = double(subs(theta_base,   physical_all, pv));
    torque_fn = cell(1, nDOF);
    for jn = 1:nDOF
        torque_fn{jn} = matlabFunction(subs(results(jn).Torque, physical_all, pv), 'Vars', {q, qd, qdd});
    end
    e_full = 0; e_max = 0; tau_max = 0;
    for k = 1:n_verification_states
        qk = -pi + 2*pi*rand(1,nDOF); qdk = -2 + 4*rand(1,nDOF); qddk = -3 + 6*rand(1,nDOF);
        tau = zeros(nDOF,1);
        for jn = 1:nDOF, tau(jn) = torque_fn{jn}(qk, qdk, qddk); end
        e_full  = max(e_full, norm(tau - A_common_fn(qk,qdk,qddk)*theta_common_num, inf));
        e_max   = max(e_max,  norm(tau - A_base_fn(qk,qdk,qddk)*theta_base_num, inf));
        tau_max = max(tau_max, norm(tau, inf));
    end
    fprintf('Verification: e_max = %.3e (max |tau| = %.3e)\n', e_max, tau_max);
    t_stage5 = toc(t_start) - t_stage12 - t_stage3 - t_stage4;

    %% ========================================================================
    % CLASSICAL COUNT: LINK DIMENSIONS SUBSTITUTED NUMERICALLY
    % ========================================================================
    % Re-express every monomial in the standard parameters
    % [XX..ZZ, MX=m*Lcx, MY=m*Lcy, MZ=m*Lcz, M=m] and fix the geometry.
    std_symbols = sym([]);  lc_symbols = sym([]);  lc_replacements = sym([]);
    for i = 1:nDOF
        mi = sym(sprintf('m%d',i), 'real');
        axes_ = 'xyz';
        for a = 1:3
            Lc  = sym(sprintf('Lc%d%s', i, axes_(a)), 'real');
            MLc = sym(sprintf('M%s%d', upper(axes_(a)), i), 'real');
            lc_symbols(end+1) = Lc; lc_replacements(end+1) = MLc/mi; %#ok<AGROW>
            std_symbols(end+1) = MLc; %#ok<AGROW>
        end
        for nm = {'Ixx','Ixy','Ixz','Iyy','Iyz','Izz'}
            std_symbols(end+1) = sym(sprintf('%s%d', nm{1}, i), 'real'); %#ok<AGROW>
        end
        std_symbols(end+1) = mi; %#ok<AGROW>
    end
    std_symbols = [std_symbols, actuation_symbols];
    theta_std = subs(theta_common, lc_symbols, lc_replacements);
    theta_std = simplify(subs(theta_std, sym(geometry_names), geometry_values));
    J_std = jacobian(theta_std, std_symbols);
    is_linear_in_standard = isempty(symvar(J_std));
    if is_linear_in_standard
        Phi_std = Phi_str * double(J_std);
        sv_std  = svd(Phi_std);
        r_std   = sum(sv_std > rank_factor*max(size(Phi_std))*max(sv_std));
    else
        r_std = NaN;   % centroidal parameterization: not linear in standard parameters
    end
    fprintf('Base count with numerical geometry: %d\n', r_std);

    %% ========================================================================
    % SAVE
    % ========================================================================
    summary = struct('robot', robot_name, 'p_i', p_i, 'p', p, 'M', M, 'r', r, ...
        'p_b', p_b, 'rank_tolerance', rank_tol, 'sigma_r', sv_str(r), ...
        'sigma_r_plus_1', sv_str(min(r+1,numel(sv_str))), 'e_full', e_full, ...
        'e_max', e_max, 'tau_max', tau_max, 'n_verification_states', n_verification_states, ...
        'r_numeric_geometry', r_std, 'run_time_s', toc(t_start), ...
        'stage_times_s', [t_stage12 t_stage3 t_stage4 t_stage5]);

    save(fullfile(results_folder, 'base_regressor.mat'), 'A_matrix_base', 'theta_base', ...
        'A_matrix_common', 'theta_common', 'independent_cols', 'dependent_cols', ...
        'Gamma', 'summary');

    fid = fopen(fullfile(results_folder, 'theta_base.txt'), 'w');
    fprintf(fid, 'theta_base = %s;\n', char(theta_base));
    fclose(fid);
    fid = fopen(fullfile(results_folder, 'theta_base_list.txt'), 'w');
    for k = 1:p_b, fprintf(fid, 'Theta_b%-3d = %s\n', k, char(theta_base(k))); end
    fclose(fid);
    fid = fopen(fullfile(results_folder, 'A_matrix_base.txt'), 'w');
    fprintf(fid, 'A_matrix_base = %s;\n', char(A_matrix_base));
    fclose(fid);

    fid = fopen(fullfile(results_folder, 'regrouping_summary.txt'), 'w');
    fprintf(fid, 'Robot                                  : %s\n', robot_name);
    fprintf(fid, 'Link dimensions                        : %s\n', ternary(opt.known_geometry, 'numeric (known geometry)', 'symbolic'));
    fprintf(fid, 'Actuation terms (Ia, Fv, Fc)           : %s\n', ternary(opt.actuation, 'included', 'not included'));
    fprintf(fid, 'Monomials per joint                    : %s\n', mat2str(p_i));
    fprintf(fid, 'Distinct monomials p                   : %d\n', p);
    fprintf(fid, 'Rank-test states M                     : %d\n', M);
    fprintf(fid, 'Structural rank r = p_b                : %d\n', p_b);
    fprintf(fid, 'Dependent directions removed           : %d (%.2f %%)\n', p-p_b, 100*(p-p_b)/p);
    fprintf(fid, 'sigma_r / sigma_(r+1) / tolerance      : %.4e / %.4e / %.4e\n', ...
        summary.sigma_r, summary.sigma_r_plus_1, rank_tol);
    fprintf(fid, 'Verification states                    : %d\n', n_verification_states);
    fprintf(fid, 'max |tau - A*Theta|_inf                : %.3e\n', e_full);
    fprintf(fid, 'max |tau - A_b*Theta_b|_inf (e_max)    : %.3e\n', e_max);
    fprintf(fid, 'max |tau|_inf                          : %.3e\n', tau_max);
    fprintf(fid, 'Base count with numerical geometry     : %d\n', r_std);
    fprintf(fid, 'Run time                               : %.1f s\n', summary.run_time_s);
    fprintf(fid, 'Stage times 1-2 / 3 / 4 / 5 (s)        : %.1f / %.1f / %.1f / %.1f\n', summary.stage_times_s);
    fclose(fid);
    type(fullfile(results_folder, 'regrouping_summary.txt'));

end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end
