function stats = misstdoa_brief_visualization(sol)
% MISSTDOA_BRIEF_VISUALIZATION Plot status report for current solution

opt.threshold = 0.15; % TODO: Options
opt.stdout = 2; % TODO: Options

[zcalc,zok] = misstdoa_calc_z(sol);
zerr = zcalc(:)-sol.z(:);
zinl = sol.inlmatrix(:);
zerr = zerr(find(zinl));

stats = [...
length(sol.rows) size(sol.z,1); ...
length(sol.cols) size(sol.z,2); ...
sum(sol.inlmatrix(:)) prod(size(sol.z)); ...
sum(zok(:)) prod(size(sol.z)); ...
norm(zerr) 0];
DETMISS=1;
DETIN=2;
ZOK=4;
STDOUT=8;
OUTLIER=16;
LAST=32;
dist = abs(zcalc-sol.z);
detin = sol.inlmatrix;
detmiss = isnan(sol.z);
plot_im = ...
    (detmiss)*DETMISS + ...       % missing data (black)
    (detin)*DETIN + ...         % true pos (green) - good
    zok*ZOK + ...
    (dist > std(zerr)*opt.stdout)*STDOUT + ...
    (dist > opt.threshold)*OUTLIER + ...
    0; % Last line
%    (~zok & detin)*3 + ... % included outliers (yellow) quite bad
%    (~detmiss & ~detin)*3;  % + ... % false pos (red) really bad

%    (truein & ~detin)*2 + ... % false neg (orange) quite bad
%    (trueout & ~detin)*4; % true neg (yellow) - good
plot_colmap = [0 1 0;0.2 0.7 0.2; 1 0 0;0.9 1 0; 0 0 0]; % G DG R Y K
plot_colmap = ones(LAST,3);
plot_colmap(DETMISS,:) = [0 0 0]; % Black
plot_colmap(DETIN,:) = [0 1 0]; % Green
plot_colmap(DETIN+ZOK,:) = [0 1 0]; % Green
plot_colmap(DETIN+STDOUT,:) = [0.2 0.7 0.2]; % Dark green
plot_colmap(DETIN+ZOK+STDOUT,:) = [0.2 0.7 0.2]; % Dark green
plot_colmap(STDOUT,:) = [0.9 1 0]; % Yellow
plot_colmap(DETIN+OUTLIER,:) = [1 0.5 0]; % Orange
plot_colmap(DETIN+STDOUT+OUTLIER,:) = [1 0.5 0]; % Orange
plot_colmap(DETIN+ZOK+OUTLIER,:) = [1 0.5 0]; % Orange
plot_colmap(DETIN+ZOK+STDOUT+OUTLIER,:) = [1 0.5 0]; % Orange
plot_colmap(STDOUT+OUTLIER,:) = [1 0 0]; % Red
plot_colmap(OUTLIER,:) = [1 0 0]; % Red
names = []

%
figure(1); clf; subplot(2,1,1);
plot_im(plot_im < 1) = LAST - plot_im(plot_im < 1);
image(plot_im);
colormap(plot_colmap);
title(['R:' num2str(length(sol.rows)) '/' num2str(size(sol.z,1)) ...
    ' C:' num2str(length(sol.cols)) '/' num2str(size(sol.z,2)) ...
    ' I:' num2str(sum(sol.inlmatrix(:))) '/' num2str(prod(size(sol.z))) ...
    ' P:' num2str(sum(zok(:))) '/' num2str(prod(size(sol.z)))  ...
    ' E: ' num2str(norm(zerr)) ' std ' num2str(std(zerr))]);
subplot(2,1,2);
h = histogram(zerr,100);

subplot(2,1,2);
plot_im = plot_im(zinl);

[m] = unique(plot_im);
for n = 1:size(m,1)
    k = m(n);
    counts(n,:) = histcounts(zerr(plot_im==k), h.BinEdges);
end
h = bar(h.BinEdges(1:end-1),counts,'stacked','FaceColor','flat', 'BarWidth', 1);
colormap(plot_colmap);
for n = 1:size(m,1)
    k = m(n);
    if k < 1 || k > LAST, k = LAST; end
    h(n).CData = plot_colmap(k,:);
end
m = length(sol.rows);
n = length(sol.cols);
title(['Theta2 - E: ' num2str(norm(zerr)*length(zerr)/(length(zerr)-(3*m+3*n-6))) ...
    ' Std: ' num2str(std(zerr)*length(zerr)/(length(zerr)-(3*m+3*n-6))) ...
        ' Avg ', num2str(mean(zerr))]);
