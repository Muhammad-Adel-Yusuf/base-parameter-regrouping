function [r, sv, Phi] = mdh_numeric_rank(mdh, actuation, n_states, seed)
%MDH_NUMERIC_RANK Base-parameter count with known geometry from the Newton-Euler reference alone.
%
%   [r, sv, Phi] = mdh_numeric_rank(mdh)
%   [r, sv, Phi] = mdh_numeric_rank(mdh, actuation, n_states, seed)
%
%   The recursive Newton-Euler reference mdh_rnea_reference is linear in the
%   standard parameters [XX XY XZ YY YZ ZZ MX MY MZ M] of every link (and in Ia,
%   Fv, Fc if actuation is true). Its regressor is built column by column from
%   unit parameter vectors at n_states random states, and the numerical rank is
%   the base-parameter count for the numerical link dimensions of the table.
%   This does not use the symbolic Lagrange torques, so it is an independent
%   check of the count that regroup_base_parameters reports as 'Base count with
%   numerical geometry'.
%
%   r    numerical rank; sv the singular values; Phi the (n*n_states)-by-p regressor.

    if nargin < 2 || isempty(actuation), actuation = false; end
    if nargin < 3 || isempty(n_states), n_states = max(60, 12*mdh.n); end
    if nargin < 4 || isempty(seed), seed = 31; end
    n = mdh.n;
    ref = mdh_rnea_reference(mdh);
    dn = mdh.dim_names;  dv = mdh.dim_values;

    rng(seed, 'twister');
    Q   = -pi + 2*pi*rand(n_states, n);
    QD  = -2 + 4*rand(n_states, n);
    QDD = -3 + 6*rand(n_states, n);

    % parameter columns: each is a list of (name, value) pairs that realises one unit standard parameter
    cols = {};
    ax = 'xyz';  ten = {'Ixx','Ixy','Ixz','Iyy','Iyz','Izz'};
    for i = 1:n
        for k = 1:6, cols{end+1} = struct('on', {{sprintf('%s%d', ten{k}, i)}}, 'val', 1, 'off', {{}}); end %#ok<AGROW>
        for a = 1:3        % first moment: (m = 1, Lc = 1) minus (m = 1, Lc = 0)
            cols{end+1} = struct('on', {{sprintf('m%d', i), sprintf('Lc%d%s', i, ax(a))}}, 'val', [1 1], ...
                'off', {{sprintf('m%d', i)}}); %#ok<AGROW>
        end
        cols{end+1} = struct('on', {{sprintf('m%d', i)}}, 'val', 1, 'off', {{}}); %#ok<AGROW>
    end
    if actuation
        for i = 1:n
            for nm = {'Ia','Fv','Fc'}
                cols{end+1} = struct('on', {{sprintf('%s%d', nm{1}, i)}}, 'val', 1, 'off', {{}}); %#ok<AGROW>
            end
        end
    end
    p = numel(cols);

    Phi = zeros(n*n_states, p);
    for s = 1:n_states
        rows = n*(s-1) + (1:n);
        for c = 1:p
            C = cols{c};
            t1 = ref(Q(s,:), QD(s,:), QDD(s,:), [dn, C.on], [dv, C.val]);
            if isempty(C.off)
                Phi(rows, c) = t1;
            else
                t0 = ref(Q(s,:), QD(s,:), QDD(s,:), [dn, C.off], [dv, 1]);
                Phi(rows, c) = t1 - t0;                 % the first-moment column is the difference
            end
        end
    end
    sv = svd(Phi);
    r = sum(sv > 1e-8*max(size(Phi))*max(sv));
end
