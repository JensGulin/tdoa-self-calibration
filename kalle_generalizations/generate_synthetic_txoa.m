function [z, gt] = generate_synthetic_txoa(m, n, dims, type, sigma, miss_ratio, out_ratio, out_range)
% GENERATE_SYNTHETIC_TDOA Generate synthetic TDOA measurements
%   [z, gt] = GENERATE_SYNTHETIC_TDOA(m, n, dim) generates an m by n matrix
%       z of TDOA measurements between m receivers and n senders. The
%       measurements satisfy z_ij = || r_i - s_j || + o_j, where r_i and
%       s_j are receiver and senders positions embedded in a
%       dim-dimensional space. gt is a struct containing ground truth data
%       for z, r, s, o, and d = z-o.
%   [z, gt] = generate_synthetic_tdoa(m, n, dim, sigma, miss_ratio, out_ratio, out_range)
%       Gaussian noise, missing data, and outliers are added to the
%       measurements. The outliers are uniformly sampled from the interval
%       out_range.

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

if length(dims)==1;
    r_dim = dims(1);
    s_dim = dims(1);
else
    r_dim = dims(1);
    s_dim = dims(2);
end
max_dim = max(r_dim,s_dim);

r = [randn(r_dim, m);zeros(max_dim-r_dim, m)];
s = [randn(s_dim, n);zeros(max_dim-s_dim, n)];
switch type
    case 'tdoa'
        o = randn(1, n);
    case 'toa'
        o = zeros(1, n);
    case 'cotoa'
        o = ones(1, n)*randn(1,1);
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
