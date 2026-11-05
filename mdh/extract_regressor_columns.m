function [theta, A, info] = extract_regressor_columns(tau, kin, dim_symbols, dim_values, known, log_fn)
%EXTRACT_REGRESSOR_COLUMNS Stage 1 on nested torque expressions (derivative route).
%
%   [theta, A, info] = extract_regressor_columns(tau, kin, dim_symbols, dim_values, known, log_fn)
%
%   tau          n-by-1 symbolic torques (physical symbols and kinematic symbols)
%   kin          1-by-3n kinematic symbols [q qd qdd]
%   dim_symbols  link-dimension symbols; dim_values their numerical values
%   known        true: substitute the dimensions before extraction
%   log_fn       optional function handle log_fn(fmt, ...)
%
%   Every torque is a polynomial in the physical symbols, so the regressor
%   column of a monomial mu = prod x_j^alpha_j is
%       c(q,qd,qdd) = (1/alpha!) d^alpha tau / dx^alpha  evaluated at x = 0.
%   The monomials are found from the coefficients of the expanded torque at two
%   random numerical kinematic states (a coefficient that vanishes at both is a
%   measure-zero event); columns that are identically zero are dropped. Each
%   joint is checked: tau - A*theta must vanish at a further random state.
%
%   theta{j}: monomials of joint j; A{j}: their columns (1-by-p_j symbolic).
%
%   One setting, sequential. extract_regressor_columns_multi does several
%   settings at once and can use a parallel pool.

    if nargin < 6 || isempty(log_fn), log_fn = @(varargin) []; end
    [th, AA, inf_] = extract_regressor_columns_multi(tau, kin, dim_symbols, dim_values, known, log_fn, false);
    theta = th{1};  A = AA{1};  info = inf_{1};
end
