function mdh = mdh_table(name, rows, dims)
%MDH_TABLE Validated modified-DH table of a serial manipulator (revolute and prismatic joints).
%
%   mdh = mdh_table(name, rows, dims)
%
%   name  robot name (used for the results folder)
%   rows  n-by-4 numeric or symbolic matrix, one row per joint i = 1..n:
%             [ alpha_{i-1}   a_{i-1}   d_i   theta0_i ]
%         in the modified Denavit-Hartenberg (Khalil-Kleinfinger / Craig)
%         convention, where the joint angle is theta_i = q_i + theta0_i and
%             {i-1}T{i} = Rot_x(alpha_{i-1}) Trans_x(a_{i-1}) Rot_z(theta_i) Trans_z(d_i).
%         Entries may contain the link-dimension symbols listed in dims
%         (for example L1, -L4, pi/2). Classic (distal) DH tables must be
%         converted first.
%   dims  struct whose field names are the link-dimension symbols used in
%         rows and whose values are the known numerical dimensions in metres,
%         for example struct('L1', 0.333, 'L2', 0.316). Use struct() if the
%         table is purely numeric.
%
%   The inertial parameters (masses, centres of mass, inertia tensors) are not
%   part of the table: they stay symbolic throughout the regrouping.
%
%   A fifth column (0 revolute, 1 prismatic) declares the joint type; without
%   it every joint is revolute. For a revolute joint the joint variable q_i
%   is added to the angle, theta_i = q_i + theta0_i, and d_i is a constant. For
%   a prismatic joint the variable is added to the offset, d_i -> q_i + d_i,
%   and theta0_i is the constant joint angle.
%
%   Returns a struct with the fields name, n, alpha, a, d, theta0 (n-by-1
%   symbolic), prismatic (n-by-1 logical), dim_names, dim_values and dim_symbols.

    if nargin < 3, dims = struct(); end
    if ~ischar(name) && ~isstring(name), error('mdh_table:name', 'name must be text.'); end
    rows = sym(rows);
    if size(rows, 2) ~= 4 && size(rows, 2) ~= 5
        error('mdh_table:columns', ['rows must have 4 columns [alpha a d theta0] (all joints revolute) ' ...
            'or 5 columns with a joint type (0 revolute, 1 prismatic); found %d.'], size(rows, 2));
    end
    n = size(rows, 1);
    if size(rows, 2) == 5
        prismatic = rows(:, 5);
        if ~all(arrayfun(@(e) isequal(e, sym(0)) || isequal(e, sym(1)), prismatic))
            error('mdh_table:type', 'Column 5 (joint type) must be 0 (revolute) or 1 (prismatic).');
        end
        prismatic = logical(double(prismatic));
        rows = rows(:, 1:4);
    else
        prismatic = false(n, 1);
    end
    if n < 2, error('mdh_table:joints', 'At least two joints are required.'); end

    dim_names = fieldnames(dims).';
    dim_values = zeros(1, numel(dim_names));
    dim_symbols = sym(zeros(1, numel(dim_names)));
    for k = 1:numel(dim_names)
        v = dims.(dim_names{k});
        if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v)
            error('mdh_table:dims', 'Dimension %s must be a finite numeric scalar.', dim_names{k});
        end
        dim_values(k) = v;
        dim_symbols(k) = sym(dim_names{k}, 'real');
    end

    % every free symbol of the table must be a declared dimension
    used = symvar(rows);
    declared = arrayfun(@char, dim_symbols, 'UniformOutput', false);
    for k = 1:numel(used)
        if ~any(strcmp(char(used(k)), declared))
            error('mdh_table:symbol', ['Symbol %s in the table is not declared in dims. ' ...
                'Add it to dims with its numerical value.'], char(used(k)));
        end
    end
    % declare the symbols as real so that they match the torque expressions
    if ~isempty(used)
        rows = subs(rows, used, arrayfun(@(k) sym(char(k), 'real'), used));
    end
    % a wrongly entered classic-DH table usually has a nonzero theta0 and
    % alpha that is not a multiple of pi/2; only the second is checked
    al = rows(:, 1);
    for i = 1:n
        if isempty(symvar(al(i)))
            k = double(al(i)) / (pi/2);
            if abs(k - round(k)) > 1e-9
                warning('mdh_table:alpha', ['Row %d: alpha = %g rad is not a multiple of pi/2. ' ...
                    'Check that the table is in the modified (proximal) convention.'], i, double(al(i)));
            end
        end
    end

    mdh = struct('name', char(name), 'n', n, 'alpha', rows(:, 1), 'a', rows(:, 2), ...
        'd', rows(:, 3), 'theta0', rows(:, 4), 'prismatic', prismatic, 'dim_names', {dim_names}, ...
        'dim_values', dim_values, 'dim_symbols', dim_symbols);
end
