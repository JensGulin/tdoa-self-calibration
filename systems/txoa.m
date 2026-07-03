function [r, s, o, sol] = txoa(z, varargin)
% TXOA Perform TXOA self-calibration using robust methods
%   [r, s, o, sol] = TXOA(z) performs TXOA self-calibration in 3D using the
%       TXOA measurements in z (m x n). The estimated receiver (3 x m) and
%       sender (3 x n) positions, and the offsets (1 x n) are returned
%       together with a solution structure containing additional
%       information. The method works as follows:
%           1. Initialize a relaxed solution in u, v, a, b, o.
%           2. Extend solution to more rows and columns.
%           3. Upgrade to solution in r, s, o.
%   [r, s, o, sol] = TXOA(z, ...) additional inputs:
%       display - controls the amount of printing.
%           off, none - no printing.
%           iter - printouts for each iteration of RANSAC/optimization.
%       sigma - estimated standard deviation of the measurement noise. This
%           is used to set suitable thresholds for RANSAC.
%       offsetsolver - a structure with the offset solver to use.
%
% See M. Larsson et al. (2021) Fast and Robust Stratified Self-Calibration
% Using Time-Difference-of-Arrival Measurements for details.
%
% See also get_offset_solvers.


%keyboard;

% Use (9r/5s) offset solver by default.
default_solver = get_offset_solvers(3, 9, 5);

% Parse inputs.
p = inputParser;
valid_display = @(x) ismember(x, {'off', 'none', 'iter'});
addParameter(p, 'display', 'off', valid_display);
%addParameter(p, 'offsetsolver', default_solver);
addParameter(p, 'offset_type', 'tdoa');
addParameter(p, 'dims', [3 3]);
addParameter(p, 'sigma', 0.01);
addParameter(p, 'mul1', 4); % Sets threshold based on sigma.
addParameter(p, 'mul2', 6); % Sets threshold based on sigma.
addParameter(p, 'mul3', 10); % Sets threshold based on sigma.
addParameter(p, 'mul4', 12); % Sets threshold based on sigma.
addParameter(p, 'mul5', inf); % Sets threshold based on sigma.
addParameter(p, 'iters', 1000);
parse(p, varargin{:});
opts = p.Results;

solver_key = sprintf('%s_rank%d', lower(opts.offset_type), min(opts.dims));
switch solver_key
    case 'tdoa_rank1'
        asolver = get_offset_solver_func(@solver_tdoa_rank1_53);
    case 'tdoa_rank2'
        asolver = get_offset_solver_func(@solver_tdoa_rank2_74);
    case 'tdoa_rank3'
        asolver = get_offset_solver_func(@solver_tdoa_rank3_95);

    case 'cotoa_rank1'
        asolver = get_offset_solver_func(@solver_cotoa_rank1_33);
    case 'cotoa_rank2'
        asolver = get_offset_solver_func(@solver_cotoa_rank2_44);
    case 'cotoa_rank3'
        asolver = get_offset_solver_func(@solver_cotoa_rank3_55);

    case 'toa_rank1'
        asolver = get_offset_solver_func(@solver_toa_rank1_32);
    case 'toa_rank2'
        asolver = get_offset_solver_func(@solver_toa_rank2_43);
    case 'toa_rank3'
        asolver = get_offset_solver_func(@solver_toa_rank3_54);

    otherwise
        error('Unknown solver case: %s', solver_key);
end
addParameter(p, 'offsetsolver', asolver);
parse(p, varargin{:});
opts = p.Results;


% TODO: Generalize for other dimensions than 3D.

%keyboard;

% Find initial u, v, a, b, o.
sol = init_uvabo_ransac(z, ...
    'display', opts.display, 'threshold', opts.mul1*opts.sigma,...
    'solver', opts.offsetsolver, 'offset_type', opts.offset_type, ...
    'rank', min(opts.dims), 'iters', opts.iters);

% Bundle over u, v, a, b, o.
sol = refine_uvabo(sol, 'display', opts.display);

% Extend solutions to more rows and columns.
sol = extend_uvabo_ransac(sol, ...
    'display', opts.display, 'threshold', opts.mul2*opts.sigma);

% Upgrade to solution in r, s, o.
sol = upgrade_ransac(sol, ...
    'display', opts.display, 'threshold', opts.mul3*opts.sigma, ...
    'dims', opts.dims);
% sol = upgrade_linear(sol);

% Extend solutions to more rows and columns.
% TODO: Clean up and rename functions below.
sol = misstdoa_reestimate_cols_rso(sol, 'threshold', opts.mul4*opts.sigma);
sol = misstdoa_reestimate_rows_rso(sol, 'threshold', opts.mul4*opts.sigma);

% Local optimization over L, q.
% sol = refine_Lq(sol);

% Local optimization over r, s, o.
sol = refine_rso(sol, 'display', opts.display, 'max_iters', 200);
if isfinite(opts.mul5)
    sol = refine_rso_robust(sol, 'display', opts.display, 'threshold', opts.mul5*opts.sigma);
end

% Output.
[m, n] = size(z);
dim = 3;
r = nan(dim, m);
s = nan(dim, n);
o = nan(1, n);

r(:, sol.rows) = sol.r;
s(:, sol.cols) = sol.s;
o(:, sol.cols) = sol.o;
