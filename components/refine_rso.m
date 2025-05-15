function sol = refine_rso(sol, varargin)
% REFINE_RSO Perform local optimization of receiver and senders
%   sol = REFINE_RSO(sol) performs local optimization over the receiver and
%       sender positions, and the offsets.
%   sol = REFINE_RSO(sol, ...) additional inputs:
%       display - controls the amount of printing.
%           off, none - no printing
%           iter - printouts for each iteration of RANSAC/optimization
%       max_iters - maximum number of interations in optimization.
%       tol - tolerance for RMS error.
% TODO JAG Merge with refine_rso_robust

% Parse inputs.
p = inputParser;
valid_display = @(x) ismember(x, {'off', 'none', 'iter'});
addParameter(p, 'display', 'off', valid_display);
addParameter(p, 'max_iters', 100);
addParameter(p, 'tol', 1e-6);
parse(p, varargin{:});
opts = p.Results;

% Display.
if ~any(strcmpi(opts.display, {'off', 'none'}))
    fprintf('Refining solution over R, S, o.\n');
end

%keyboard;

% Extract variables for convenience.
z = sol.z(sol.rows, sol.cols);
r = sol.r;
s = sol.s;
o = sol.o;
inliers = sol.inlmatrix(sol.rows, sol.cols);

[I, J] = find(inliers);
Z = z(inliers);

for i = 1:opts.max_iters
    [res, jac] = calcresandjac(r, s, o, I, J, Z);

    %%
    m1 = size(r,2);
    n1 = size(s,2);
    N = 3 * m1 + 3 * n1 + n1;  % TODO JAG MAGIC dimensions?
    Jo = (1:n1) + N - n1;  % TODO JAG name and explain
    dontmoveindex = []; 

    % TODO JAG reindent and check. tdoa?
        EE = speye(N, N);
        % First fix the offset indices
        switch sol.offset_type
            case 'cotoa'
                EE(Jo,Jo(1)) = ones(n1,1); % This is so that all offsets are changed in the same way
                dontmoveindex = [dontmoveindex Jo(2:end)]; % Remove the other columns.
            case 'toa'
                dontmoveindex = [dontmoveindex Jo(1:end)]; % Remove the other columns.
        end
        EE(:, dontmoveindex) = [];

    if strcmpi(opts.display, 'iter')
        fprintf('Iter %3d: norm(res)=%e, rms(res)=%e\n', i, norm(res), rms(res));
    end

    if rms(res) < opts.tol
        break;
    end

    % Gauss-Newton step.
    jac0 = jac * EE;
    dz = -EE*((jac0' * jac0 + 1e-4 * speye(size(jac0, 2))) \ (jac0' * res)); % TODO MAGIC. See tol?

    [rnew, snew, onew] = update(r, s, o, dz);
    resnew = calcresandjac(rnew, snew, onew, I, J, Z);

    [norm(res) norm(res+jac*dz) norm(resnew)]  % TODO JAG conditional or discard?

    % If no improvement, try reducing the step size.
    j = 0;
    while norm(resnew) >= norm(res)
        dz = dz / 2;
        [rnew, snew, onew] = update(r, s, o, dz);
        resnew = calcresandjac(rnew, snew, onew, I, J, Z);
        j = j + 1;
        if j > 50 % TODO MAGIC
            break;
        end
    end

    if norm(resnew) < norm(res)
        r = rnew;
        s = snew;
        o = onew;
    else
        if strcmpi(opts.display, 'iter')
            fprintf('Stalled\n');
        end
        break;
    end
end

sol.r = r;
sol.s = s;
sol.o = o;
end

function [res, jac] = calcresandjac(r, s, o, I, J, Z)
v = r(:, I) - s(:, J);
d = vecnorm(v);
res = d' + o(J)' - Z;

if nargout > 1
    dim = size(r, 1);
    m = size(r, 2);
    n = size(s, 2);
    nres = length(Z);

    vid = v ./ d;
    II = repelem((1:nres)', 2*dim+1, 1);
    JJ = [dim * I + (1 - dim:0), dim * m + dim * J + (1 - dim:0), dim * (m + n) + J]';
    VV = [vid; -vid; ones(1, nres)];
    jac = sparse(II(:), JJ(:), VV(:), nres, dim*(m + n)+n);
end
end

function [r, s, o] = update(r, s, o, dz)
dim = size(r, 1);
m = size(r, 2);
n = size(s, 2);

dzr = dz(1:dim*m);
dzs = dz(dim*m+1:dim*(m + n));
dzo = dz(dim*(m + n)+1:dim*(m + n)+n);

r(:) = r(:) + dzr;
s(:) = s(:) + dzs;
o(:) = o(:) + dzo;
end
