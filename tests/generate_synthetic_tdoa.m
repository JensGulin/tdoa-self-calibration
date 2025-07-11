function [z, gt] = generate_synthetic_tdoa(m, n, dims, varargin)
% GENERATE_SYNTHETIC_TDOA Generate synthetic TDOA measurements
%   [z, gt] = GENERATE_SYNTHETIC_TDOA(m, n, dims) generates an m by n matrix
%       z of TxOA measurements between m receivers and n senders. The
%       measurements satisfy z_ij = || r_i - s_j || + o_j, where r_i and
%       s_j are receiver and senders positions embedded in a 'dims'-dimensional space. 
%       'dims' is the number of dimensions (2 or 3) for positions or
%       a two-number array [r_dim s_dim], if different.
%       gt is a struct containing ground truth data for z, r, s, o, and d = z-o.
%   [z, gt] = GENERATE_SYNTHETIC_TDOA(m, n, dims, sigma, miss_ratio, out_ratio, out_range)
%       Gaussian noise, missing data, and outliers are added to the
%       measurements. The outliers are uniformly sampled from the interval
%       out_range. Noise is present in both gt.z and in output z.
%       The gt.gt_z field carries the accurate z in this case.
%  See <a href="matlab:help GENERATE_SYNTHETIC_TXOA">generate_synthetic_txoa</a>, this is simply a wrapper:
%       TxOA type (here 'tdoa') can then be 'toa' or 'cotoa' instead,
%       and controls the offset o.

[z, gt] = generate_synthetic_txoa(m, n, dims, 'tdoa', varargin{:});
