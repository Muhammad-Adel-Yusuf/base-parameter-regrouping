function torques = load_torque_equations(folder, nDOF)
%LOAD_TORQUE_EQUATIONS Read Torque_1.txt ... Torque_n.txt as symbolic expressions.
%   Each file contains MATLAB assignments ending with "Torque_k = <expr>;".
%   Files may also be split into parts (Ak_1 = ...; Torque_k = Ak_1 + ...;).
%   Every identifier in the files is declared as a real symbol first.

    if nargin < 2, nDOF = 6; end
    torques = cell(1, nDOF);
    excluded = {'sin','cos','tan','sqrt','exp','log','abs','sign','pi', ...
                'matrix','Inf','NaN'};
    for k = 1:nDOF
        fname = fullfile(folder, sprintf('Torque_%d.txt', k));
        if ~isfile(fname), error('Cannot find %s', fname); end
        txt = fileread(fname);
        % char() of a 1x1 symbolic array may print as matrix([[...]])
        txt = regexprep(txt, 'matrix\(\[\[\s*', '(');
        txt = regexprep(txt, '\s*\]\]\)', ')');
        torques{k} = eval_torque_text(txt, sprintf('Torque_%d', k), excluded);
    end
end

function T = eval_torque_text(txt, target, excluded)
    ids = unique(regexp(txt, '\<[A-Za-z][A-Za-z0-9_]*\>', 'match'));
    ids = setdiff(ids, excluded);
    for i = 1:numel(ids)
        if isvarname(ids{i})
            eval(sprintf('%s = sym(''%s'', ''real'');', ids{i}, ids{i}));
        end
    end
    eval(txt);
    T = eval(target);
end
