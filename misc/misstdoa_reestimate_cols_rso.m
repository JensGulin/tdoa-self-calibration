function [solut] = misstdoa_reestimate_cols_rso(sol, varargin)
% [solut]=misstdoa_reestimate_cols_theta1(sol);
% TODO JAG merge with misstdoa_reestimate_rows_rso?

% Parse inputs.
p = inputParser;
valid_display = @(x) ismember(x, {'off', 'none', 'iter'});
addParameter(p, 'display', 'off', valid_display);
addParameter(p, 'iters', 10);
addParameter(p, 'threshold', 0.1);
parse(p, varargin{:});
opts = p.Results;

% Display.
if ~any(strcmpi(opts.display, {'off', 'none'}))
    fprintf('Reestimating sender positions and offsets\n');
end

z = sol.z;
r = sol.r;
rows = sol.rows;

for jj = 1:size(z,2)
    %% Check early that column has enough data for trilateration
    % TODO JAG MAGIC <5 <4 nr_inliers
    switch sol.offset_type
        case 'tdoa'
            if sum(isfinite(z(rows, jj)))<5, % Changed to 5 got error in trilateration
                continue;
            end
            [yy, oo, inlny, nr_inliers] = tdoa_trilateration_y_one_ransac(z(rows, jj), r, opts.iters, opts.threshold);
            if nr_inliers < 4
                continue;
            end
            yy = real(yy);
            oo = real(oo);
            [sny, ony, resny] = tdoa_trilateration_y_one_bundle(z(rows, jj), r, yy, oo, inlny);
            nr_inliers_ny = length(inlny);
            inlny2 = false(length(rows), 1);  % TODO JAG flip to align with the others?
            inlny2(inlny) = true;
            rms_ny = sqrt(resny'*resny);
        case 'cotoa'
            oo = median(sol.o);
            if sum(isfinite(z(rows, jj)))<4, % Changed to 4 got error in trilateration
                continue;
            end
            [yy, inlny, nr_inliers] = toa_trilateration_one_ransac(z(rows,jj)-oo, r, opts.iters, opts.threshold);
            if nr_inliers < 3
                continue;
            end
           %[yy,oo,inlny,nr_inlierid,err_rms] = tdoa_trilateration_y_one_ransac(z(rows,jj),r,ransac_k,ransac_tol);
           yy = real(yy);
           [sny, resny] = toa_trilateration_one_bundle(z(rows,jj)-oo, r, yy, inlny);
           ony = oo;
           nr_inliers_ny = length(inlny);
           inlny2 = false(1, length(rows));
           inlny2(inlny) = true;
           rms_ny = sqrt(resny'*resny);
        case 'toa'
            oo = 0;
            if sum(isfinite(z(rows, jj)))<4, % Changed to 4 got error in trilateration
                continue;
            end
            [yy, inlny, nr_inliers] = toa_trilateration_one_ransac(z(rows,jj), r, opts.iters, opts.threshold);
            if nr_inliers < 3
                continue;
            end
           %[yy,oo,inlny,nr_inlierid,err_rms] = tdoa_trilateration_y_one_ransac(z(rows,jj),r,ransac_k,ransac_tol);
           yy = real(yy);
           [sny, resny] = toa_trilateration_one_bundle(z(rows,jj), r, yy, inlny);
           ony = oo;
           nr_inliers_ny = length(inlny);
           inlny2 = false(1, length(rows));
           inlny2(inlny) = true;
           rms_ny = sqrt(resny'*resny);
    end
    %[length(inlny) std(resny)]

    % Is there a previous guess?
    tmp = find(sol.cols == jj);
    if length(tmp) < 1,
        % No previous estimate of this column
        tmp = size(sol.cols, 2) + 1;
        sol.s(:, tmp) = sny;
        sol.o(:, tmp) = ony;
        sol.cols(tmp) = jj;
        sol.inlmatrix(rows, jj) = inlny2;
    elseif length(tmp) == 1,
        % There is a previous estimate
        s0 = sol.s(:, tmp);
        o0 = sol.o(:, tmp);
        inl0 = sol.inlmatrix(:, jj);

        zproj = tdoa_calc_u_from_xyo(r, s0, o0);
        res0 = zproj - z(rows, jj);
        res0 = res0(find(sol.inlmatrix(rows, jj)));
        nr_inliers_0 = sum(inl0);
        rms_0 = sqrt(res0'*res0);
        %[sum(inl0) length(res0) std(res0)]

        if (nr_inliers_ny > nr_inliers_0) || ((nr_inliers_ny == nr_inliers_0) && (rms_ny < rms_0))
            sol.s(:, tmp) = sny;
            sol.o(:, tmp) = ony;
            sol.inlmatrix(rows, jj) = inlny2;
        end
    else
        % Det h?r vore konstigt.
        error('');
    end
end

solut = sol;


if 0,
    if nargin < 3,
        y0 = NaN * ones(x_dim, n);
    end

    if nargin < 4,
        o0 = NaN * ones(1, n);
    end

    if nargin < 5,
        index = find(ones(1, n));
    end;

    if nargin < 6,
        inliers = isfinite(u);
    end;

    inl = zeros(size(inliers));

    ransac_tol = 0.02;
end;
