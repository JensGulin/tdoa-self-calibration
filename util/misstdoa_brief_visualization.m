function stats = misstdoa_brief_visualization(sol, varargin)
% MISSTDOA_BRIEF_VISUALIZATION Plots status report for a solution
%  That is, the solution is reprojected to TxOA and compared with given z.
%  1. The first part shows the (color coded) status of each z separately.
%  The stats shows the number of processed columns (C, mics)
%  and rows (R, events). I is the number of inliers, P the processed z
%  (C*R). The error E is total squared error (Note: not mean) and std the
%  error standard deviation (not squared). The error only includes the inliers.
%  2. Second part is the histogram of inlier errors, with additional
%  measurements shown if they are in view (e.g. true inliers falsely marked
%  as outliers).
%  The stats here are the same E and std as before, but (somehow) scaled.
%  In addition the (unscaled) mean error (not squared) indicates any bias.
%  3. The color coding is the same in both parts, with legend given in the
%  second part. The distinction of a true outlier (TN) and a false outlier
%  (FN) is the 'threshold' option. The separation of true inliers (TP) and
%  questionable inliers (?P) is controlled by the 'stdin' option.
% 
%  stats = MISSTDOA_BRIEF_VISUALIZATION(sol)
%  stats = MISSTDOA_BRIEF_VISUALIZATION(sol, varargin)
% Input:
%  * sol: The solution struct (in e.g. rso or uvabo form).
%  * varargin: Optional options. It can be a mixture of 
%    positional (unnamed) parameters in the order below, a struct
%    with named fields, or ('name', value) parameters. If the same
%    parameter is given several times, the last one wins. Positional
%    parameters must come before any others and do not accept strings,
%    use string array or cell array for the positional title. 
% Options:
%  * fig: [int] The figure number. Def. 1.
%  * threshold: [double] The error accepted for a true inlier. Def. 0.15.
%  * stdin: [double] The error expected for a true inlier, in terms of
%    number of std deviations from zero. Def. 2.
%  * title: [string] Figure title. Can be several lines if string array.
%    Uses ('Interpreter', 'none'), but all properties can be changed by
%    submitting a cell array, e.g. {"title",'FontSize',10}. Def. "".
% Output:
%  * stats: An array of (5,2) size with various numbers.

% Parse inputs.
p = inputParser;
% p.KeepUnmatched = true; % No error for unknown
addOptional(p, 'fig', 1);
addOptional(p, 'threshold', 0.15);
addOptional(p, 'stdin', 2);
addOptional(p, 'title', "");
parse(p, varargin{:});
opts = p.Results;

if ~isfield(sol,'rows') || ~isfield(sol,'cols')
    [sol.rows, sol.cols] = size(sol.z);
    sol.type = 'z';
end

[zcalc,zok] = misstdoa_calc_z(sol);
zerr_all = zcalc(:)-sol.z(:);
zinl = sol.inlmatrix(:);
zerr = zerr_all(zinl);

stats = [...
length(sol.rows) size(sol.z,1); ...
length(sol.cols) size(sol.z,2); ...
sum(sol.inlmatrix(:)) numel(sol.z); ...
sum(zok(:)) numel(sol.z); ...
norm(zerr) 0];
DETMISS=1;
DETIN=2;
NOK=4;
STDOUT=8;
OUTLIER=16;
LAST=32;
LEFTOUT=STDOUT+DETMISS;
dist = abs(zcalc-sol.z);
detin = sol.inlmatrix;
detmiss = isnan(sol.z);
plot_im = ...
    (detmiss)*DETMISS + ...       % missing data (black)
    (detin)*DETIN + ...         % true pos (green) - good
    ~zok*NOK + ...  % not evaluated
    (dist > std(zerr)*opts.stdin)*STDOUT + ...
    (dist > opts.threshold)*OUTLIER + ...
    0; % Last line
plot_im(plot_im == 0) = (dist(plot_im == 0) <= std(zerr)*opts.stdin)*(LEFTOUT);

plot_colmap = [0 1 0;0.2 0.7 0.2; 1 0 0;0.9 1 0; 0 0 0]; % G DG R Y K
plot_colmap = ones(LAST,3); % White
plot_colmap(DETMISS,:) = [0 0 0]; % Black
plot_colmap(DETMISS + NOK,:) = [0 0 0]; % Black
plot_colmap(NOK,:) = [0.5 0.5 0.5]; % Gray
plot_colmap(DETIN,:) = [0 1 0]; % Green
plot_colmap(DETIN+STDOUT,:) = [0.2 0.7 0.2]; % Dark green
plot_colmap(STDOUT,:) = [0.9 1 0]; % Yellow
plot_colmap(NOK+STDOUT,:) = [0.9 1 0]; % Yellow
plot_colmap(LEFTOUT,:) = [0.9 1 0]; % Yellow
plot_colmap(DETIN+OUTLIER,:) = [1 0.5 0]; % Orange
plot_colmap(DETIN+STDOUT+OUTLIER,:) = [1 0.5 0]; % Orange
plot_colmap(DETIN+NOK,:) = [1 0.5 0]; % Orange
plot_colmap(DETIN+NOK+STDOUT,:) = [1 0.5 0]; % Orange
plot_colmap(DETIN+NOK+OUTLIER,:) = [1 0.5 0]; % Orange
plot_colmap(DETIN+NOK+STDOUT+OUTLIER,:) = [1 0.5 0]; % Orange
plot_colmap(NOK+OUTLIER,:) = [1 0 0]; % Red
plot_colmap(NOK+STDOUT+OUTLIER,:) = [1 0 0]; % Red
plot_colmap(STDOUT+OUTLIER,:) = [1 0 0]; % Red
plot_colmap(OUTLIER,:) = [1 0 0]; % Red
plot_colmap(LEFTOUT,:) = [0.93,0.93,0.7]; % [0 0 1]; % Blue


plot_colname = ones(LAST,1); % White
plot_colname(DETMISS,:) = 7; % Black
plot_colname(DETMISS + NOK,:) = 7; % Black
plot_colname(NOK,:) = 8; % Gray
plot_colname(DETIN,:) = 2; % Green
plot_colname(DETIN+STDOUT,:) = 3; % Dark green
plot_colname(STDOUT,:) = 4; % Yellow
plot_colname(NOK+STDOUT,:) = 4; % Yellow
plot_colname(LEFTOUT,:) = 4; % Yellow
plot_colname(DETIN+OUTLIER,:) = 5; % Orange
plot_colname(DETIN+STDOUT+OUTLIER,:) = 5; % Orange
plot_colname(DETIN+NOK,:) = 5; % Orange
plot_colname(DETIN+NOK+STDOUT,:) = 5; % Orange
plot_colname(DETIN+NOK+OUTLIER,:) = 5; % Orange
plot_colname(DETIN+NOK+STDOUT+OUTLIER,:) = 5; % Orange
plot_colname(NOK+OUTLIER,:) = 6; % Red
plot_colname(NOK+STDOUT+OUTLIER,:) = 6; % Red
plot_colname(STDOUT+OUTLIER,:) = 6; % Red
plot_colname(OUTLIER,:) = 6; % Red


%
figure(opts.fig); clf;
if (isstring(opts.title) || ischar(opts.title)) && all(opts.title ~= "")
    sgtitle(opts.title, 'Interpreter', 'none');
elseif iscell(opts.title)
    sgtitle(opts.title{:});
end
subplot(2,1,1);
plot_im(plot_im < 1) = LAST - plot_im(plot_im < 1);
image(plot_im);
colormap(plot_colmap);
title(['R:' num2str(length(sol.rows)) '/' num2str(size(sol.z,1)) ...
    ' C:' num2str(length(sol.cols)) '/' num2str(size(sol.z,2)) ...
    ' I:' num2str(sum(sol.inlmatrix(:))) '/' num2str(numel(sol.z)) ...
    ' P:' num2str(sum(zok(:))) '/' num2str(numel(sol.z))  ...
    ' E: ' num2str(norm(zerr)) ' std ' num2str(std(zerr))]);

%%
subplot(2,1,2);
% Make a first plot to get h.BinEdges
h = histogram(zerr_all,1000);
%h = histogram(zerr,100);

cname = {"white: unknown", "green: inlier (ok, TP)", "dark: inlier (tail, ?P)", "yellow: inlier (excluded, FN)","orange: bad (included, FP)","red: outlier (excluded, TN)","black: missing (excluded)", "gray: unused (excluded, ??)"};
%plot_im = plot_im(zinl); % TODO: Do we want to keep only inliers?

labels = {};
[m] = unique(plot_im);
counts = zeros(size(m,1),size(h.BinEdges,2)-1);
for n = 1:size(m,1)
    k = m(n);
    counts(n,:) = histcounts(zerr_all(plot_im==k), h.BinEdges);

end
h = bar(h.BinEdges(1:end-1),counts,'stacked','FaceColor','flat', 'BarWidth', 1);
colormap(plot_colmap);
for n = 1:size(m,1)
    k = m(n);
    if k < 1 || k > LAST, k = LAST; end
    h(n).CData = plot_colmap(k,:);
    labels(end+1) = {num2str(k) + " " + cname{plot_colname(k)}};
end
legend(labels);
m = length(sol.rows);
n = length(sol.cols);
% TODO JAG: What is Theta2 and the scaling
% *length(zerr)/(length(zerr)-(3*m+3*n-6))) ??
% In paper https://ieeexplore-ieee-org.ludwig.lub.lu.se/stamp/stamp.jsp?tp=&arnumber=9414309
% Theta2 is (page 3) the UVabo estimation, but I don't see the relation,
% unless it's left-overs from a time when the function only did that?
title(['Theta2 - E: ' num2str(norm(zerr)*length(zerr)/(length(zerr)-(3*m+3*n-6))) ...
    ' Std: ' num2str(std(zerr)*length(zerr)/(length(zerr)-(3*m+3*n-6))) ...
        ' Avg ', num2str(mean(zerr))]);
