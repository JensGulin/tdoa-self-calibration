function solout = upgrade_ransac(sol, varargin)
% UPGRADE_RANSAC Upgrade a relaxed solution
%   sol = UPGRADE_RANSAC(sol) upgrades a relaxed solution, a solution in u,
%       v, a, b and o, to a solution in r, s and o. This is done using
%       robust methods (see note below).
%   sol = UPGRADE_RANSAC(sol, ...) additional inputs:
%       display - controls the amount of printing.
%           off, none - no printing.
%           iter - printouts for each iteration of RANSAC/optimization.
%       iters - number of iterations in RANSAC loop.
%       threshold - threshold for the absolute error in TDOA measurment
%           when classifying inliers/outliers.
%
% See M. Larsson et al. (2020) Upgrade Methods for Stratified Sensor
% Network Self-calibration for details.

% Parse inputs.
p = inputParser;
valid_display = @(x) ismember(x, {'off', 'none', 'iter'});
addParameter(p, 'display', 'off', valid_display);
addParameter(p, 'iters', 50);
addParameter(p, 'threshold', 0.1);
addParameter(p, 'dims', [3 3]);
parse(p, varargin{:});
opts = p.Results;

if length(opts.dims)==1,
    opts.dims = opts.dims*ones(1,2);
end

% Display.
if ~any(strcmpi(opts.display, {'off', 'none'}))
    fprintf('Upgrading solution using RANSAC.\n');
end

% Initialize output struct.
solout.type = 'rso';
solout.rows = sol.rows;
solout.cols = sol.cols;
solout.inlmatrix = sol.inlmatrix;
solout.z = sol.z;
solout.o = sol.o;

% Convert sol to the problem structure used in the previous paper.
prob = struct();
prob.a = sol.b;
prob.b = sol.a;
prob.c = 0;
prob.U = sol.u';
prob.V = sol.v;
prob.Dmeas = sol.z(sol.rows, sol.cols) - sol.o;
inliers = sol.inlmatrix(sol.rows, sol.cols);

% Options controlling local optimization in minimal solvers.
solvopts = struct();
solvopts.tol = 1e-9;
solvopts.maxIters = 20;
solvopts.refine = false;

% Get all upgrade solvers.
dimtype = [num2str(opts.dims(1)) num2str(opts.dims(2))];

switch dimtype
    case '33'
        solvers = getSolvers();
    case '22'
        solvers = getSolvers_2D2D();
    case '23'
        solvers = getSolvers_2D3D();
    case '32'
        solvers = getSolvers_3D2D();
end
% Keep only those who work with the number of receivers and senders.
allids = vertcat(solvers.id);
keep = length(sol.rows) >= allids(:, 1) + 1 & length(sol.cols) >= allids(:, 2) + 1;
solvers = solvers(keep);

%keyboard;

% RANSAC loop.
max_inliers = 0;
for i = 1:opts.iters
    % Pick a solver randomly.
    % TODO: Is this the best strategy?
    solver = solvers(randi(length(solvers)));

    % Create random minimal sample from the problem. Since u and v are
    % dense we do not need to worry about missing data.
    sample = createRandomSample(solver.id, prob);

    switch dimtype
        case '32'
            % transpose sample
            sample_transpose = sample;
            sample_transpose.U = sample.V;
            sample_transpose.V = sample.U;
            sample_transpose.a = sample.b';
            sample_transpose.b = sample.a';
            sample_transpose.indU = sample.indV;
            sample_transpose.indV = sample.indU;
            sample_transpose.fullU = sample.fullV;
            sample_transpose.fullV = sample.fullU;
            % solve minimal problem
            [Lti, q] = solver.solve(sample_transpose, solvopts);
        otherwise
            % Solve minimal problem.
            [Lti, q] = solver.solve(sample, solvopts);
    end


    for j = 1:length(Lti)
        switch dimtype
            case '33'
                Rhat = Lti{j} * sample.fullU;
                Shat = Lti{j}' \ (sample.fullV + q(:, j));
            case '22'
                Rhat = Lti{j} * sample.fullU;
                Shat = Lti{j}' \ (sample.fullV + q(:, j));
            case '23'
                Rhat = Lti{j} * sample.fullU;
                Shat = Lti{j}' \ (sample.fullV + q(:, j));
                Dhat = pdist2(Rhat', Shat');
                S3_squared = nanmedian(prob.Dmeas.^2-Dhat.^2);
                Rhat = [Rhat;zeros(1,size(Rhat,2))];
                Shat = [Shat;sqrt(relu(S3_squared))];
            case '32'
                RhatT = Lti{j} * sample_transpose.fullU;
                ShatT = Lti{j}' \ (sample_transpose.fullV + q(:, j));
                % transpose solution
                Rhat = ShatT;
                Shat = RhatT;
                % Add third dimension
                Dhat = pdist2(Rhat', Shat');
                R3_squared = nanmedian(prob.Dmeas.^2-Dhat.^2,2)';
                Rhat = [Rhat;sqrt(relu(R3_squared))];
                Shat = [Shat;zeros(1,size(Shat,2))];
        end

        % Calculate error in distances.
        Dhat = pdist2(Rhat', Shat');
        err = abs(prob.Dmeas(inliers)-Dhat(inliers));
        ninl = sum(err < opts.threshold);

        if ninl > max_inliers
            max_inliers = ninl;
            solout.r = Rhat;
            solout.s = Shat;
            solout.L = inv(Lti{j}');
            solout.q = q(:, j);
            solout.u = sample.fullU';
            solout.v = sample.fullV;
            solout.offset_type = sol.offset_type;
            solout.dims = opts.dims;

            if strcmpi(opts.display, 'iter')
                fprintf('Iter %3d: inliers = %3d, solver = %s\n',...
                    i, ninl, solver.name);
            end
        end
    end
end

if max_inliers == 0
    error('RANSAC upgrade failed to find a solution.');
end
