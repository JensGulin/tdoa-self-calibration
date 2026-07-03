function [bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, varargin)
% INIT_UVABO_RANSAC Initialize a relaxed TDOA solution
%   sol = INIT_UVABO_RANSAC(z) initializes a relaxed solution, i.e., a
%       solution in u, v, a, b and o, using robust methods.
%   sol = INIT_UVABO_RANSAC(z) additional inputs:
%           display - controls the amount of printing.
%               off, none - no printing.
%               iter - printouts for each iteration of RANSAC/optimization.
%           iters - number of iterations in RANSAC loop.
%           solver - an array of offset solvers to use.
%           threshold - threshold for the absolute error in TDOA measurment
%               when classifying inliers/outliers.

% 1. The offsets are solved for first using the minimal offset solvers.
% 2. The minimal solution is then extended to include more columns using
%    robust methods.
% 3. Inliers are counted and the best solution is kept.

%keyboard;

% Default solver.
default_solver.solv = @solver_tdoa_rank3_95;
default_solver.name = func2str(default_solver.solv);
default_solver.m = 9;
default_solver.n = 5;
default_solver.rank = 3;

% Parse inputs.
p = inputParser;
valid_display = @(x) ismember(x, {'off', 'none', 'iter'});
addParameter(p, 'display', 'off', valid_display);
addParameter(p, 'iters', 100000);
addParameter(p, 'solver', default_solver);
%addParameter(p, 'threshold', 0.2);
addParameter(p, 'threshold', 0.01);
addParameter(p, 'offset_type', 'tdoa');
addParameter(p, 'rank', 3);
parse(p, varargin{:});
opts = p.Results;

% %TODO: Reconsider solver selection, read rank etc from the solver.
% % If no solver given, use the function to get them.
% % Tried to do this now
% solver_key = sprintf('%s_rank%d', lower(opts.offset_type), opts.rank);
% switch solver_key
%     case 'tdoa_rank1'
%         asolver = get_offset_solver_func(@solver_tdoa_rank1_53);
%     case 'tdoa_rank2'
%         asolver = get_offset_solver_func(@solver_tdoa_rank2_74);
%     case 'tdoa_rank3'
%         asolver = get_offset_solver_func(@solver_tdoa_rank3_95);
% 
%     case 'cotoa_rank1'
%         asolver = get_offset_solver_func(@solver_cotoa_rank1_33);
%     case 'cotoa_rank2'
%         asolver = get_offset_solver_func(@solver_cotoa_rank2_44);
%     case 'cotoa_rank3'
%         asolver = get_offset_solver_func(@solver_cotoa_rank3_55);
% 
%     case 'toa_rank1'
%         asolver = get_offset_solver_func(@solver_toa_rank1_32);
%     case 'toa_rank2'
%         asolver = get_offset_solver_func(@solver_toa_rank2_43);
%     case 'toa_rank3'
%         asolver = get_offset_solver_func(@solver_toa_rank3_54);
% 
%     otherwise
%         error('Unknown solver case: %s', solver_key);
% end
% addParameter(p, 'solver', asolver);
% parse(p, varargin{:});
% opts = p.Results;

%keyboard;


% Display.
if ~any(strcmpi(opts.display, {'off', 'none'}))
    fprintf('Running RANSAC over U, V, a, b, o.\n');
end

% Init RANSAC.
bestsol = [];
max_inliers = 0;
best_err = Inf;
stats1 = zeros(1, opts.iters);
stats2 = zeros(3, 0);
stats2_counter = 0;

% The problem type defines how much data is needed,
% based on how many unknowns there are.
% Typically 4-6 depending on rank and offset_type.
% Note: Additional data is needed to verify the minimal solution,
% so we usually count inliers excluding the data used by the solver.
switch opts.offset_type
    case 'tdoa'
        minsolver_nr_rows = opts.rank + 2;
    case {'toa','cotoa'}
        minsolver_nr_rows = opts.rank + 1;
end
% TODO: Reconsider this. Suppose the criteria is > minimal for
% verification and the value is dictated by the solver. If the solver is
% smaller than minsolver_nr_rows it should be discarded, right?


for iRansac = 1:opts.iters
    % m receivers and n senders.
    [m, n] = size(z);
    ok = isfinite(z);
    
    % Pick random solver.
    solver = opts.solver(randi(length(opts.solver)));

    % m_solv, n_solv are the number of rows and columns needed for the minimal
    % solver.
    m_solv = solver.m;
    n_solv = solver.n;
    rank_solv = solver.rank;

    % Precompute matrices needed for going from (z,o) to (u,v,a,b);
    wr = zeros(m_solv, 1);
    wr(1) = 1;
    ws = zeros(n_solv, 1);
    ws(1) = 1;

    Cr = eye(m_solv) - wr * ones(1, m_solv);
    Cs = eye(n_solv) - ws * ones(1, n_solv);
    
    % Choose m_solv random rows.
    rowsol = randperm(m, m_solv);

    % Select all columns which has data for all rows.
    ok_cols = find(all(ok(rowsol, :), 1));
    if length(ok_cols) < n_solv
        continue;
    end
    colsol = ok_cols(randperm(length(ok_cols), n_solv));

    sols = solver.solv(z(rowsol, colsol));
    nsols = size(sols, 2);
    stats1(1, iRansac) = nsols;
    for iSol = 1:nsols
        osol = sols(:, iSol).';
        if any(abs(imag(osol))./abs(osol) > 1e-6)
            continue; % Discard imaginary solutions
        end
        osol = real(osol);
        zsol = z(rowsol, colsol);
        dsol = zsol - osol;
        if any(dsol(:) <= 0) % TODO: Need some margin here?
            continue;
        end

        % TODO: Local optimization over offsets?

        % Find initial estimate of u,v,a,b from zsol.
        dsol2 = dsol.^2;
        M = Cr' * dsol2 * Cs / (-2); % This M is M/(-2) compared to the ICASSP 2020 paper
        [uu, ss, vv] = svd(M);
        u = uu(:, 1:rank_solv);
        v = ss(1:rank_solv, 1:rank_solv) * vv(:, 1:rank_solv)';
        a = dsol2 * ws - wr' * dsol2 * ws;
        b = wr' * dsol2;

        % The following should now hold:
        % dsol.^2 = (zsol-osol).^2 == -2*(u*v)+a+b

        % TODO: Clean up code below and extract function. There is overlap
        % with extend_vbo_ransac, but otherwise fine to keep in here?

        % Now check for inliers on the rest of the cols (on the rows of the
        % solution). Try all cols where a subset of the solution rows satisfy;
        % (i) there are strictly more than "minsolver_nr_rows" measurements, and
        okcol = find(sum(ok(rowsol, :)) > minsolver_nr_rows); % Det här kanske man ska modifiera (done)
        % (ii) col was not in solution already
        restcols = setdiff(okcol, colsol);
        inliersrest = zeros(size(restcols));
        v_rest = zeros(rank_solv, length(restcols));
        b_rest = zeros(1, length(restcols));
        o_rest = zeros(1, length(restcols));
        inl_rest = zeros(m_solv, length(restcols));
        nr_inliers = 0;
        tot_err = 0;
        for iCol = 1:length(restcols)
            %keyboard;
            col_idx = restcols(iCol);
            % For each new column,
            % use the fact that we know (u,a)
            % and calculate (vny,bny,ony)

            % First select a minimum number of ok rows for the solver.
            okrow = find(ok(rowsol, col_idx));
            rows_cut = okrow(randperm(length(okrow), minsolver_nr_rows));
            z_cut = z(rowsol(rows_cut), col_idx);
            u_cut = u(rows_cut, :);
            a_cut = a(rows_cut);

            % Solve the unknowns from,
            % -2*(u*[v])+a+[b] == (z-[o]).^2
            switch opts.offset_type
                case 'tdoa'
                    % In this case we need to get ony as well.
                    % Expand (z-ony).^2 as z.^2 -2z.*ony + ony.^2
                    % and rearrange to
                    % -2u*[vny] -1*[ony.^2] + 2z*[ony] + 1*[bny] == z.^2 - a
                    AAA = [(-2*u_cut), -ones(minsolver_nr_rows, 1),... 
                        2*z_cut, ones(minsolver_nr_rows, 1)];
                    bbb = z_cut.^2 - a_cut;
                    x_part = AAA \ bbb;
                    ony = x_part(end-1);
                    % Adjust xxx if o^2 was off the mark.
                    lamb = ony^2 - x_part(end-2);
                    x_hom = [zeros(1,opts.rank), 1, 0, 1]';
                    xxx = x_part + lamb * x_hom;
                case 'cotoa'
                    ony = osol(1);
                    d2_cut = (z_cut-ony).^2;
                    AAA = [(-2*u_cut),ones(minsolver_nr_rows, 1)];
                    bbb = d2_cut - a_cut;
                    xxx = AAA \ bbb;
                case 'toa'
                    % (z(cc,col_idx)-ony).^2  is equal to -2*(u*vny)+a+bny
                    ony = 0;
                    d2_cut = (z_cut-ony).^2;
                    AAA = [(-2*u_cut),ones(minsolver_nr_rows, 1)];
                    bbb = d2_cut - a_cut;
                    xxx = AAA \ bbb;
            end
            vny = xxx(1:rank_solv);
            bny = xxx(end);
            v_rest(:, iCol) = vny;
            b_rest(1, iCol) = bny;
            o_rest(1, iCol) = ony;
            % ... and then check for inliers among measurements in
            % this column
            err = (sqrt(relu(-2 * (u(okrow, :) * vny) + a(okrow) + bny)) + ony) - z(rowsol(okrow), col_idx);
            inlid = find(abs(err) < opts.threshold);
            if length(inlid) > minsolver_nr_rows
                inliersrest(iCol) = 1;
                inl_rest(okrow(inlid), iCol) = ones(length(inlid), 1);
                nr_inliers = nr_inliers + length(inlid) - minsolver_nr_rows;
                tot_err = tot_err + sum(err(inlid).^2); % I am adding the five zeros here, but nevermind.
            end
        end


        tot_err = sqrt(tot_err);
%         nr_inliers = sum(inliersrest) + n_solv;
        nr_inliers = nnz(inl_rest) + m_solv * n_solv;

        stats2_counter = stats2_counter + 1;
        stats2(1, stats2_counter) = iRansac;
        stats2(2, stats2_counter) = nr_inliers;
        stats2(3, stats2_counter) = tot_err;
        stats2(4,stats2_counter)=osol(1);

        if (nr_inliers > max_inliers) || ((nr_inliers == max_inliers) && (tot_err < best_err))
            if strcmpi(opts.display, 'iter')
                fprintf('Iter %5d: inliers = %3d, error = %e\n', iRansac, nr_inliers, tot_err);
            end
            %keyboard;

            max_inliers = nr_inliers;
            best_err = tot_err;

            inliersrest = find(inliersrest);
            bestsol.rows = rowsol;
            bestsol.cols = [colsol, restcols(inliersrest)];
            bestsol.row1 = bestsol.rows(1);
            bestsol.col1 = bestsol.cols(1);

            %keyboard;
            o = [osol, o_rest(inliersrest)];
            v = [v, v_rest(:, inliersrest)];
            b = [b, b_rest(inliersrest)];

            bestsol.inlmatrix = false(m, n);
            bestsol.inlmatrix(bestsol.rows, bestsol.cols) = [true(m_solv, n_solv), inl_rest(:, inliersrest)];
            bestsol.z = z;
            bestsol.o = o;
            bestsol.b = b;
            bestsol.a = a;
            bestsol.u = u;
            bestsol.v = v;
            bestsol.type = 'uvabo';
            bestsol.offset_type = opts.offset_type;
            bestsol.rank = opts.rank;
        end
    end
end

if isempty(bestsol)
    error('Failed to find a solution.');
end

%keyboard;

if strcmpi(opts.display, 'iter') % TODO: Option to enable plots?
    misstdoa_briefer_report(bestsol);
    misstdoa_brief_visualization(bestsol);
end
