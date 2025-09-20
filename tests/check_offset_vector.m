function [ok, err] = check_offset_vector(sols_o, gt_o, varargin)
% CHECK_OFFSET_VECTOR Verify that there is an acceptable o-vector match.
% [ok, err] = CHECK_OFFSET_VECTOR(sols_o,gt_o,opts)
%
% [ok, err] = CHECK_OFFSET_VECTOR(sol,gt,opts)
%   Detailed explanation goes here

% Parse inputs.
p = inputParser;
valid_display = @(x) ismember(x, {'off', 'none', 'iter'});
addParameter(p, 'display', 'iter', valid_display);
addParameter(p, 'tol', 5e-4);
parse(p, varargin{:});
opts = p.Results;

print = opts.display == "iter";
if isstruct(sols_o)
    % Input is the solution struct, get o
    sol = sols_o;
    sols_o = nan(size(sol.z,2),1);
    sols_o(sol.cols) = sol.o;
end
if isstruct(gt_o)
    % Input is the gt struct, get o
    gt = gt_o;
    gt_o = nan(size(gt.z,2),1);
    gt_o(gt.cols) = gt.o;
end

% Accept any transpose of the o vector
assert(min(size(gt_o)) == 1, "gt_o must be a vector or a gt struct.")
gt_o = reshape(gt_o,[],1);
assert(size(sols_o,1) == size(gt_o,1), ...
    "sols_o and gt_o need matching row count (%i vs %i).", size(sols_o,1), size(gt_o,1))

% Find cols with acceptable error
% Note: abs takes imaginary part too, which is accepted if small enough for tolerance...
err = max(abs(sols_o-gt_o));
ok = (err) < opts.tol;

if print
    disp("check_offset_vector [sols_o,gt_o]:");
    disp([sols_o,gt_o]);
    disp("Worst error (per column):");
    disp(err);
    err = min(err);
    disp("Error with precision (and threshold):");
    fprintf("    %e (%e)\n", err, opts.tol);
end

ok = nnz(ok); % Number of fully acceptable solutions
if(ok == 0)
    error("No solutions matching GT!");
end
err = sols_o-gt_o;
