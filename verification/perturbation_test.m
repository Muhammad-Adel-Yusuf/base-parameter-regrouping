%% perturbation_test.m
% Sensitivity of the numerical check (Stage 5, 200 random states) and of the
% exact certificate to an error in one dependency coefficient. The certified
% coefficients Gamma_ex are perturbed in one entry by a relative amount delta;
% the numerical residual max|A_D*Theta_D - A_I*Gamma_p*Theta_D| (random
% parameter values, as in Stage 5) is compared with the exact verdict
% (polynomial reduction, as in exact_verification.m).
% Writes ../perturbation_test.txt.
here = fileparts(mfilename('fullpath'));  root = fileparts(here);
names = {'UR5e', 'IRB120'};      % the two six-joint robots, as reported in the paper
cases = cellfun(@(r) fullfile(root, r, 'results', 'known_geometry'), names, 'UniformOutput', false);
deltas = [1e-6 1e-9 1e-12 1e-14];
fid = fopen(fullfile(root, 'perturbation_test.txt'), 'w');
fprintf(fid, '%-7s %-10s %-22s %-22s %s\n', 'robot', 'delta', 'numerical residual', 'unperturbed residual', 'exact verdict');
for c = 1:numel(cases)
    S = load(fullfile(cases{c}, 'base_regressor.mat'), 'A_matrix_common', 'independent_cols', 'dependent_cols', 'Gamma');
    A = S.A_matrix_common;  ind = S.independent_cols;  dep = S.dependent_cols;
    n = size(A, 1);
    q   = arrayfun(@(k) sym(sprintf('t%d',   k), 'real'), 1:n);
    qd  = arrayfun(@(k) sym(sprintf('t%dd',  k), 'real'), 1:n);
    qdd = arrayfun(@(k) sym(sprintf('t%ddd', k), 'real'), 1:n);
    Gex = sym(round(S.Gamma*1e12))/sym(10)^12;            % certified candidate
    A_fn = matlabFunction(A, 'Vars', {q, qd, qdd});
    rng(4242);
    X = cell(200, 1);
    for k = 1:200, X{k} = A_fn(-pi + 2*pi*rand(1,n), -2 + 4*rand(1,n), -3 + 6*rand(1,n)); end
    thD = 0.1 + rand(numel(dep), 1);
    [kmax, idx] = max(abs(double(Gex(:))));               % perturb the largest coefficient
    resid = @(G) max(cellfun(@(Ak) norm(Ak(:, dep)*thD - Ak(:, ind)*(G*thD), inf), X));
    r0 = resid(double(Gex));
    for d = deltas
        Gp = Gex;  Gp(idx) = Gp(idx)*(1 + sym(d));
        r = resid(double(Gp));
        nz = exact_remainders(A, ind, dep, Gp);
        verdict = 'certified';
        if nz > 0, verdict = sprintf('rejected (%d entries nonzero)', nz); end
        fprintf(fid, '%-7s %-10.0e %-22.3e %-22.3e %s\n', names{c}, d, r, r0, verdict);
    end
end
fclose(fid);
type(fullfile(root, 'perturbation_test.txt'));

function bad = exact_remainders(A, ind, dep, Gr)
    n = size(A, 1);
    qn = arrayfun(@(k) sym(sprintf('t%d', k), 'real'), 1:n);
    s  = arrayfun(@(k) sym(sprintf('s%d', k), 'real'), 1:n);
    c  = arrayfun(@(k) sym(sprintf('c%d', k), 'real'), 1:n);
    E = subs(expand(A(:, dep) - A(:, ind)*Gr), [sin(qn) cos(qn)], [s c]);
    vars = [s c setdiff(symvar(E), [s c])];
    bad = 0;
    for k = 1:numel(E)
        if isequal(E(k), sym(0)), continue; end
        r = polynomialReduce(E(k), s.^2 + c.^2 - 1, vars, 'MonomialOrder', 'lexicographic');
        if ~isequal(expand(r), sym(0)), bad = bad + 1; end
    end
end
