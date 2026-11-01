function export_latex_base_list(results_folder, out_file)
%EXPORT_LATEX_BASE_LIST Write theta_base as a LaTeX align block in the
% standard-parameter notation XX_i..ZZ_i, MX_i, MY_i, MZ_i, M_i.
%
%   export_latex_base_list('results', 'results/theta_base_latex.tex')

    fp = sympref('FloatingPointOutput', false);      % print exact fractions
    restore = onCleanup(@() sympref('FloatingPointOutput', fp));
    S = load(fullfile(results_folder, 'base_regressor.mat'), 'theta_base');
    th = S.theta_base(:);
    v = symvar(th);  names = arrayfun(@char, v, 'UniformOutput', false);
    % Lc_i{x,y,z}*m_i -> MX_i, MY_i, MZ_i
    old = sym([]);  new = sym([]);
    for i = 1:6
        mi = sym(sprintf('m%d', i), 'real');
        for a = 'xyz'
            n = sprintf('Lc%d%s', i, a);
            if any(strcmp(names, n))
                old(end+1) = sym(n, 'real'); %#ok<AGROW>
                new(end+1) = sym(sprintf('M%s%d', upper(a), i), 'real') / mi; %#ok<AGROW>
            end
        end
    end
    th = expand(subs(th, old, new));
    % one inline-math paragraph per parameter: TeX breaks long expressions
    % after + and - signs, which keeps two-column layouts free of overfull lines
    fid = fopen(out_file, 'w');
    fprintf(fid, '{\\small\\setlength{\\parindent}{0pt}\\raggedright\n');
    allv = symvar(th);
    masses = allv(startsWith(arrayfun(@char, allv, 'UniformOutput', false), 'm'));
    for k = 1:numel(th)
        e = rationalize(th(k));
        cand = {simplify(e, 'Steps', 60), collect(e, masses), e};
        lens = cellfun(@(c) strlength(string(latex(c))), cand);
        [~, best] = min(lens);
        s = strrep(tidy(latex(cand{best})), '\frac', '\tfrac');
        if k == numel(th), p = '.'; else, p = ';'; end
        fprintf(fid, '\\hangindent=2.2em $\\Theta_{b,%d}=%s$%s\\par\n', k, s, p);
    end
    fprintf(fid, '}\n');
    fclose(fid);
    fprintf('wrote %s (%d parameters)\n', out_file, numel(th));
end

function e = rationalize(e)
% replace floating-point coefficients (from Gamma) by exact fractions
    e = expand(e);
    if isequal(e, sym(0)), return; end
    v = symvar(e);
    [cf, tm] = coeffs(e, v);           % numeric coefficients and pure monomials
    out = sym(0);
    for k = 1:numel(cf)
        [nn, dd] = rat(double(cf(k)), 1e-9);
        out = out + sym(nn)/sym(dd) * tm(k);
    end
    e = out;
end

function s = tidy(s)
    map = {'Ixx','XX'; 'Iyy','YY'; 'Izz','ZZ'; 'Ixy','XY'; 'Ixz','XZ'; 'Iyz','YZ'};
    for k = 1:size(map,1)
        s = regexprep(s, ['\\mathrm\{' map{k,1} '\}_\{(\d)\}'], [map{k,2} '_{$1}']);
        s = regexprep(s, ['\\mathrm\{' map{k,1} '(\d)\}'], [map{k,2} '_{$1}']);
    end
    s = regexprep(s, '\\mathrm\{M([XYZ])\}_\{(\d)\}', 'M$1_{$2}');
    s = regexprep(s, '\\mathrm\{M([XYZ])(\d)\}', 'M$1_{$2}');
    s = regexprep(s, '\\mathrm\{LL\}_\{(\d)\}', 'L_{\\mathrm{e}$1}');
    s = regexprep(s, '\\mathrm\{LL(\d)\}', 'L_{\\mathrm{e}$1}');
    s = regexprep(s, '(?<![A-Za-z\\])m_\{(\d)\}', 'M_{$1}');
    s = regexprep(s, '\\,', '');
    s = strrep(s, '\left(', '(');  s = strrep(s, '\right)', ')');
end
