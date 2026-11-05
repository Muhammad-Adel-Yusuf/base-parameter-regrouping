function fh = regressor_function(expr, vars, block_chars)
%REGRESSOR_FUNCTION Numerical function of a large symbolic regressor.
%
%   fh = regressor_function(expr, vars)
%   fh = regressor_function(expr, vars, block_chars)
%
%   Same as matlabFunction(expr, 'Vars', vars), but for the large regressor
%   of the Panda. Generating one function for the whole matrix needs more
%   than 50 GB of memory, so the code is generated row by row in blocks of
%   nonzero entries of at most block_chars characters (default 2e5), each
%   written to a temporary file without the common-subexpression
%   optimization. Constant entries are stored as numbers.
%
%   fh(q, qd, qdd) with 1 x n rows returns the n x p matrix. With N x n
%   inputs (one state per row) it returns the n x p x N array of all N
%   regressors in one call, which avoids the per-call overhead of the large
%   generated code.

    if nargin < 3, block_chars = 2e5; end
    [n, p] = size(expr);
    folder = tempname;  mkdir(folder);  addpath(folder);
    blocks = struct('row', {}, 'cols', {}, 'fn', {});
    constants = zeros(0, 3);                     % [row col value]
    for i = 1:n
        cols = find(arrayfun(@(e) ~isequal(e, sym(0)), expr(i,:)));
        is_const = arrayfun(@(e) isempty(symvar(e)), expr(i,cols));
        for c = cols(is_const)
            constants(end+1,:) = [i c double(expr(i,c))]; %#ok<AGROW>
        end
        cols = cols(~is_const);
        len = arrayfun(@(e) strlength(string(char(e))), expr(i,cols));
        k = 1;
        while k <= numel(cols)
            last = k;
            while last < numel(cols) && sum(len(k:last+1)) <= block_chars, last = last + 1; end
            [~, name] = fileparts([tempname '.m']);
            fn = matlabFunction(expr(i, cols(k:last)), 'Vars', vars, ...
                'File', fullfile(folder, name), 'Optimize', false);
            blocks(end+1) = struct('row', i, 'cols', cols(k:last), 'fn', fn); %#ok<AGROW>
            k = last + 1;
        end
    end
    fh = @(varargin) evaluate_blocks(blocks, constants, n, p, varargin{:});
end

function A = evaluate_blocks(blocks, constants, n, p, varargin)
    N = size(varargin{1}, 1);
    A = zeros(n, p, N);
    for b = 1:numel(blocks)
        V = blocks(b).fn(varargin{:});           % N x k (an entry may depend on few variables)
        k = numel(blocks(b).cols);
        if size(V, 1) ~= N, V = V + zeros(N, k); end
        A(blocks(b).row, blocks(b).cols, :) = reshape(V.', 1, k, N);
    end
    for c = 1:size(constants, 1)
        A(constants(c,1), constants(c,2), :) = constants(c,3);
    end
end
