function check_solution(sol, gt, varargin)
% CHECK_SOLUTION Plots and checks a solution against GT.
%  That is, the solution is aligned to GT and error is compared with 
%  given tolerance. If requirement isn't met, and error is thrown.
%
%  MISSTDOA_BRIEF_VISUALIZATION(sol,gt)
%  MISSTDOA_BRIEF_VISUALIZATION(sol, NaN)
%  MISSTDOA_BRIEF_VISUALIZATION(sol, NaN, varargin)
% Input:
%  * sol: The solution struct (in rso form).
%  * gt: The GT struct, input NaN to disable.
%  * varargin: Optional options. It can be a mixture of 
%    positional (unnamed) parameters in the order below, a struct
%    with named fields, or ('name', value) parameters. If the same
%    parameter is given several times, the last one wins. Positional
%    parameters must come before any others and do not accept strings,
%    use cell array for the positional 'register'. 
% Options:
%  * display: [string] Control the log output. Silent unless "iter". Def. "iter".
%  * tol: [double] The error accepted, either one number or
%    a vector for [r,s,o] error separately. Inf or NaN disables
%    that particular check. Def. 5e-4.
%  * fig: [int] The figure number. Def. 2.
%  * register: [string] Which points to register (align) to. Can be 
%    "r","s","rs" or "". Def. "rs".
% Output:
%  * None.

%% Parse inputs.
p = inputParser;
valid_display = @(x) ismember(x, {'off', 'none', 'iter'});
addOptional(p, 'display', 'iter', valid_display);
addOptional(p, 'tol', 5e-4);
addOptional(p, 'fig', 2);
addOptional(p, 'register', "rs");
parse(p, varargin{:});
opts = p.Results;
% Accept NaN for 'gt'.
if isnumeric(gt) && isnan(gt)
    gt = struct('r',nan(3,1),'s',nan(3,1),'o',nan);
    opts.register = "";
    opts.tol = NaN;
end
% Accept {'rs'} for 'register'.
% TODO: Seems to work for display, is validation helping?
if iscell(opts.register)
    opts.register = opts.register{:};
end
% Expand 'tol' if needed.
opts.tol(end+1:3) = opts.tol(end);
assert(sol.type == "rso");

%% Ensure consistent matrices (pad with zeros to get 3D)
% Also sort columns and rows back
% TODO: Missing rows are zeroed, should be NaN?
% How about the extra dim?
[m,n] = size(sol.z);
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

%% Register (align = rotate + translate) result against ground truth
switch opts.register
    case "r" % (mics only)
        [Q, t] = kabsch([r], [gt.r]);
    case "rs" % (mics and sources)
        [Q, t] = kabsch([r, s], [gt.r, gt.s]);
    case "s" % (sources only)
        [Q, t] = kabsch([s], [gt.s]);
    case "" % (leave as is)
        Q = eye(size(s,1));
        t = zeros(size(s,1),1);
    otherwise % what do you want?
        error("Unknown registration: " + opts.register);
end
r = Q * r + t;
s = Q * s + t;

%% Plot results.
if opts.fig > 0
    figure(opts.fig); clf;
    plot_data = gt.r;
    plot3(plot_data(1,:),plot_data(2,:),plot_data(3,:),'bo', ...
        'DisplayName',"Mic. GT");
    hold on
    plot_data = gt.s;
    plot3(plot_data(1,:),plot_data(2,:),plot_data(3,:),'mo', ...
        'DisplayName',"Src. GT");
    plot_data = r;
    plot3(plot_data(1,:),plot_data(2,:),plot_data(3,:),'bx', 'DisplayName',"Mic. Est.")
    plot_data = plot_data + 0.1;
    text(plot_data(1,:),plot_data(2,:),plot_data(3,:), ...
        string(1:size(plot_data,2)), 'FontSize', 10, 'Color', 'blue');  % Annotate mics
    plot_data = s;
    plot3(plot_data(1,:),plot_data(2,:),plot_data(3,:),'m-x', 'DisplayName',"Src. Est.");
    plot_data = plot_data + 0.05;
    text(plot_data(1,:),plot_data(2,:),plot_data(3,:), ...
        string(1:size(plot_data,2)), 'FontSize', 7, 'Color', 'magenta');  % Annotate mics
    hold off
    axis equal
    legend();
    title("Position estimates in 3D")
end

%% Calculate results.
er = vecnorm(r-gt.r);
es = vecnorm(s-gt.s);
eo = abs(o-gt.o);
err = [rms(er,"omitnan"), rms(es,"omitnan"), rms(eo,"omitnan");
       rms(er), rms(es), rms(eo)];

if opts.display == "iter"
    fprintf('RMSE - r: %e, s: %e, o: %e\n', err');
    fprintf('MAXE - r: %e, s: %e, o: %e\n', max(er), max(es), max(eo));
    fprintf('tolE - r: %e, s: %e, o: %e\n', opts.tol);
end
%% Check
% TODO JAG: Also check NaN in r but not gt
oo = 'rso';
ooo = err(1,:) > opts.tol;
if any(ooo)
    error("check_solution: Too large RMS! Failed for: " + oo(ooo));
end
