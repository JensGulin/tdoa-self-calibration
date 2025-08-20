function check_solution(sol, gt, varargin)
%% Parse inputs.
p = inputParser;
valid_display = @(x) ismember(x, {'off', 'none', 'iter'});
addParameter(p, 'display', 'iter', valid_display);
addParameter(p, 'tol', 5e-4);
addParameter(p, 'fig', 2);
parse(p, varargin{:});
opts = p.Results;

%% Ensure consistent matrices (pad with zeros to get 3D)
% Also sort columns and rows back
% TODO: Missing rows are zeroed, should be NaN?
% How about the extra dim?
[m,n] = size(gt.z);
plotsol = sol;
oo = gt.r;
oo = [oo; zeros(abs([3 0] - size(oo)))];
gt.r = oo;
oo = gt.s;
oo = [oo; zeros(abs([3 0] - size(oo)))];
gt.s = oo;
ooo = plotsol.r;
oo = zeros(3,m);
oo(1:size(ooo,1),plotsol.rows) = ooo;
r = oo;
ooo = plotsol.s;
oo = zeros(3,n);
oo(1:size(ooo,1),plotsol.cols) = ooo;
s = oo;
o = plotsol.o;
oo = NaN(1,n);
oo(1,plotsol.cols) = o;
o = oo;

%% Register (align = rotate + translate) result against ground truth (mics only)
if 0 % TODO: options
[Q, t] = kabsch([r], [gt.r]);

r = Q * r + t;
s = Q * s + t;
end

%% Register (align = rotate + translate) result against ground truth (mics and sources)
[Q, t] = kabsch([r, s], [gt.r, gt.s]);

r = Q * r + t;
s = Q * s + t;


%% Plot results.
if opts.fig > 0
    figure(opts.fig);
    plot3(gt.r(1, :), gt.r(2, :), gt.r(3, :), 'ro');
    hold on
    plot3(gt.s(1, :), gt.s(2, :), gt.s(3, :), 'go');
    plot3(r(1, :), r(2, :), r(3, :), 'r*');
    plot3(s(1, :), s(2, :), s(3, :), 'g*');
    axis equal
    hold off
end

%% Calculate results.
er = vecnorm(r-gt.r);
es = vecnorm(s-gt.s);
eo = abs(o-gt.o);

if opts.display == "iter"
    fprintf('RMSE - r: %e, s: %e, o: %e\n', rms(er,"omitnan"), rms(es,"omitnan"), rms(eo,"omitnan"));
    fprintf('MAXE - r: %e, s: %e, o: %e\n', max(er), max(es), max(eo));
end
%% Check
if any([rms(er), rms(es), rms(eo)] > opts.tol)
    error("check_solution: Too large RMS!")
end
