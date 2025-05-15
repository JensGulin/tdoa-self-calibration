function [z, gt] = generate_synthetic_txoa(m, n, dims, type, sigma, miss_ratio, out_ratio, out_range)
% GENERATE_SYNTHETIC_DATA Generate synthetic TxOA measurements
%   [z, gt] = GENERATE_SYNTHETIC_DATA(m, n, dims) generates an m by n matrix
%       z of TxOA measurements between m receivers and n senders. The
%       measurements satisfy z_ij = || r_i - s_j || + o_j, where r_i and
%       s_j are receiver and senders positions embedded in a
%       'dims'-dimensional space. gt is a struct containing ground truth data
%       for z, r, s, o, and d = z-o.
%       'dims' is a the number of dimensions (2 or 3) for positions or
%         a two-dimensional array [r_dims s_dims] if different.
%   [z, gt] = GENERATE_SYNTHETIC_DATA(m, n, dims, type, sigma, miss_ratio, out_ratio, out_range)
%       TxOA data type (default 'tdoa') can be 'toa' or 'cotoa' instead.
%       Gaussian noise, missing data, and outliers are added to the
%       measurements. The outliers are uniformly sampled from the interval
%       out_range.
% TODO JAG merge with generate_synthetic_tdoa and change name to data?
% TODO JAG manage input

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

if isscalar(dims)
    r_dim = dims;
    s_dim = dims;
else
    r_dim = dims(1);
    s_dim = dims(2);
end
if (~isnumeric(dims) || min(r_dim,s_dim) < 2)
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
        error("Illegal type '%s'.",type)
end

d = pdist2(r', s');
z = d + o;

missing = rand(m, n) < miss_ratio;
outliers = rand(m, n) < out_ratio;

gt.r = r;
gt.s = s;
gt.o = o;
gt.d = d;
gt.z = z;
gt.inliers = ~missing & ~outliers;

z = z + sigma * randn(m, n);
z(missing) = nan;

out_low = out_range(1);
out_high = out_range(2);
z(outliers) = (out_high - out_low) * rand(nnz(outliers), 1) + out_low;
