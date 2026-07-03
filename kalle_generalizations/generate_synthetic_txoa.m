function [z, gt] = generate_synthetic_txoa(m, n, dims, type, sigma, miss_ratio, out_ratio, out_range, seed)
% GENERATE_SYNTHETIC_TXOA Generate synthetic TxOA measurements
%   [z, gt] = GENERATE_SYNTHETIC_TXOA(m, n, dims) generates an m by n matrix
%       z of TxOA measurements between m receivers and n senders. The
%       measurements satisfy z_ij = || r_i - s_j || + o_j, where r_i and
%       s_j are receiver and senders positions embedded in a 'dims'-dimensional space. 
%       'dims' is the number of dimensions (2 or 3) for positions or
%       a two-number array [r_dim s_dim], if different.
%       gt is a struct containing ground truth data for z, r, s, o, and d = z-o.
%   [z, gt] = GENERATE_SYNTHETIC_TXOA(m, n, dims, type, sigma, miss_ratio, out_ratio, out_range, seed)
%       Gaussian noise, missing data, and outliers are added to the
%       measurements. The outliers are uniformly sampled from the interval
%       out_range. Noise is present in both gt.z and in output z.
%       The gt.gt_z field carries the accurate z in this case.
%       TxOA type (default 'tdoa') can be 'toa' or 'cotoa' instead,
%       and controls the offset o.
%       If seed is supplied, rng(seed) is called before generating data.

if nargin < 4
    type = 'tdoa';
end
if nargin < 5
    sigma = 0;
end
if nargin < 6
    miss_ratio = 0;
end
if nargin < 7
    out_ratio = 0;
end
if nargin < 8
    out_range = [-2 6];
end
if nargin < 9
    seed = [];
end
if nargin >= 9 && ~isempty(seed)
    rng(seed);
end

if ~isscalar(sigma)
    error("Parameter sigma must be a number, perhaps call 'generate_synthetic_txoa' directly.");
end

if ~isnumeric(dims) || isempty(dims)
    error("Illegal dims: [%s] ", join(string(dims)))
end

r_dim = dims(1);
s_dim = dims(end);

if min(r_dim,s_dim) < 1
    error("Illegal dims: [%s] ", join(string(dims)))
end
max_dim = max(r_dim,s_dim);
assert(max_dim <= 3, "Oversize dims: [%s] ", join(string(dims)))

r = [randn(r_dim, m); zeros(max_dim-r_dim, m)];
s = [randn(s_dim, n); zeros(max_dim-s_dim, n)];
switch type
    case 'tdoa'
        o = randn(1, n);
    case 'toa'
        o = zeros(1, n);
    case 'cotoa'
        o = ones(1, n)*randn(1,1);
    otherwise
        error("Illegal type '%s', use 'tdoa', 'cotoa' or 'toa'.",type)
end

d = pdist2(r', s');
z = d + o;

missing = rand(m, n) < miss_ratio;
outliers = ~missing & (rand(m, n) < out_ratio);

gt.r = r;
gt.s = s;
gt.o = o;
gt.d = d;
gt.gt_z = z;
gt.inlmatrix = ~missing & ~outliers;
gt.missing = missing;
gt.outliers = outliers;
% Make the struct a 'solution' struct
gt.offset_type = type;
gt.type = 'gt_rso';
gt.dims = [r_dim s_dim];
gt.rank = min(gt.dims);
gt.rows = 1:m;
gt.cols = 1:n;

%% Add the different types of noise.
z = z + sigma * randn(m, n);
z(missing) = nan;

out_low = out_range(1);
out_high = out_range(2);
z(outliers) = (out_high - out_low) * rand(nnz(outliers), 1) + out_low;
%% Keep the noise as z
gt.z = z;
gt.sigma = sigma;
gt.miss_ratio = miss_ratio;
gt.out_ratio = out_ratio;
gt.out_range = out_range;
gt.seed = seed;
