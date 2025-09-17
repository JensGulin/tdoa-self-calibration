% Experiment: find best o for uvabo, without the costly
% "upgrade each suggested o" ransac

%% Set up the data if needed
% Choose the case you want below.
type = "cotoa";
do_permutaions = true;  % Compare all permutations of a single input (numeric stability)
do_permutaions = false; % Take different input squares from a larger one

%% CASE 3 - Generate CLEANISH synthetic data. COTOA RANK 3
m = 15;
n = 30;
dim = 3;
sigma = 1e-8;
miss_ratio = 0.00;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, type, sigma, miss_ratio, out_ratio);

%% CASE 3b - Generate NOISY synthetic data. COTOA RANK 3
m = 15;
n = 30;
dim = 3;
sigma = 1e-3; % NOTE: This is rather high
miss_ratio = 0.00;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, type, sigma, miss_ratio, out_ratio);

%% CASE 3c - Generate OUTLIER synthetic data. COTOA RANK 3
m = 15;
n = 30;
dim = 3;
sigma = 1e-3;
miss_ratio = 0.00;
out_ratio = 0.30;
out_range = [-1 1]*0.5;

[z, gt] = generate_synthetic_txoa(m, n, dim, type, sigma, miss_ratio, out_ratio, out_range);


%% CASE 4 - Generate CLEAN synthetic data. COTOA RANK 2
m = 15;
n = 30;
dim = 2;
sigma = 1e-8;
sigma = 1e-5;
miss_ratio = 0.00;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, type, sigma, miss_ratio, out_ratio);

%%

%%%%%%%
%%%%%%% Choose the combinatorics to try
%%%%%%% 

mx = dim + 2; % choose (solver square)
nx = mx;
sx = mx - 1; % solutions with this solver
oSolver = get_offset_solvers(dim,mx,nx, type);

if do_permutaions
    fx = mx; % all permutations from only the small square
    C = perms(1:mx);
else
    fx = 10; % sample from a larger square
    C = nchoosek(1:fx,mx);
end


%% Now run for the given C
% from z, assuming rank, find all o (no voting, only gathering o)
% sampling minimal [mx,nx] input from first fx ch and fx events.
do_center = true;
do_center = false;
center = 0;
disp("Started")
tic;
jag = [];
o1 = size(C,1);
o2 = o1;
jag = nan(1+mx+nx,o2*o1*sx);
for j = 1:o2
for n=1:o1
    if any(isnan(z(C(n,:),C(j,:))),'all')
        continue; 
    end
    zx = z(C(n,:),C(j,:));
    if do_center
        center = mean(zx,'all');
        zx = zx - center;
    end
    [allsols] = oSolver.solv(zx);
    if size(allsols,2) == 0
        continue; 
    end
    allsols = allsols + center;
    jag(1,(sx*n+sx*o1*(j-1)-sx) + (1:sx)) = allsols(1,:);
    jag(2:end,(sx*n+sx*o1*(j-1)-sx) + (1:sx)) = repmat([C(n,:),C(j,:)]',1,sx);
end % n
end % j
t = toc;
fprintf("Done after %f s\n", t)
size(jag)
jag(2+mx:end,:) = jag(2+mx:end,:) + 50;
allsols = jag;

jag = allsols(:,imag(allsols(1,:)) == 0);
size(jag)
jag = jag(:,abs(jag(1,:)) < 120);
%oo = sort(,2);
size(jag)
fprintf("Took %f s (%i/%i) avg %1e ms\n", t, j,o1, 1000*t/j);

%% Show first histogram
figure(216);
%clf
h = histogram(jag(1,:),o1*o2*sx);
[a,b] = sort(h.Values);
j1 = a(end-5:end);
j2 = b(end-5:end);
[h.BinEdges(j2)],[j1]
j1 = a(end);
j2 = b(end);
[h.BinEdges(j2:j2+1)], [j1, size(jag)]

if 0
j1 = jag(1,:);
j1 = (j1 > h.BinEdges(j2)) & (j1 < h.BinEdges(j2+1));
h = histogram(jag(1,j1),o1*o2*sx);
[j1,j2] = max(h.Values);
[h.BinEdges(j2:j2+1)], [j1, size(jag)]
end 

%% TODO: Just a test
%oo = jag(:,abs(jag(1,:) - 2) < .5);
%figure(315);clf;
%h = histogram(oo(1,:),o1*o2*sx);
%layered_bar_plot(oo, 2:11, 30, 316);
%oo = jag(:,abs(jag(1,:) - 3) < .5);
%layered_bar_plot(oo, 2:11, 30, 315);

%% Decide which peak to study
j2 = floor(mean(b(end-1:end)))
j2 = b(end-1)
j2 = b(end-3)
j2 = b(end)
h.Values(j2)
j1 = h.BinEdges(j2);
j2 = h.BinEdges(j2+1);

%% Reasonable o, based on actual z
fprintf("max zx: %.20f, max z: %.20f\n", max(z(1:fx,1:fx),[],'all'),max(z,[],'all'));
fprintf("min zx: %.20f, min z: %.20f\n", min(z(1:fx,1:fx),[],'all'),min(z,[],'all'));
fprintf("avg zx: %.20f, avg z: %.20f\n", mean(z(1:fx,1:fx),'all','omitnan'),mean(z,'all','omitnan'));
fprintf("med zx: %.20f, med z: %.20f\n", median(z(1:fx,1:fx),'all','omitnan'),median(z,'all','omitnan'));

fprintf("avg zx err: %.20e, avg z err: %.20e\n", mean(z(1:fx,1:fx) - gt.gt_z(1:fx,1:fx),'all','omitnan'),mean(z - gt.gt_z,'all','omitnan'));
fprintf("sd zx err: %.20e, sd z err: %.20e\n", std(z(1:fx,1:fx) - gt.gt_z(1:fx,1:fx),0,'all','omitnan'),std(z - gt.gt_z,0,'all','omitnan'));
%% Peak
fprintf("hist peak: %.20f - %.20f (gt %.20f diff %.20e )\n",j1,j2,gt.o(1),gt.o(1)-j1)
nnn = mean(jag(1,:));
fprintf("avg all: %.20e (gt %.20f diff %.20e )\n",nnn,gt.o(1),gt.o(1)-nnn)
nnn = median(jag(1,:));
fprintf("med all: %.20e (gt %.20f diff %.20e )\n",nnn,gt.o(1),gt.o(1)-nnn)
oo = jag(:,(jag(1,:) >= j1) & ((jag(1,:)) <= j2));
assert(size(oo,2) > 0);
nnn = mean(oo(1,:));
fprintf("avg peak: %.20e (gt %.20f diff %.20e )\n",nnn,gt.o(1),gt.o(1)-nnn)
nnn = size(oo,1);
[bc,be] = layered_bar_plot(oo, 2:nnn, 30, 1316);
%% Are there any z(m,n)-reading that is under-represented (i.e outlier) for this o?
[j1,j2]
oo = jag(:,(jag(1,:) >= j1) & ((jag(1,:)) <= j2));
good = zeros(o1,o2);
for i = 1:size(oo,2)
    good(oo(2:mx+1,i), oo(mx+2:end,i)-50) = good(oo(2:mx+1,i), oo(mx+2:end,i)-50) +1;
end
figure(315);clf;
imagesc(size(oo,2) - good(1:fx,1:fx));
title("Peak count of bad? hits")
colorbar();

%% The details are however more nuanced
[~,j2] = max(bc)
j1 = be(j2);
j2 = be(j2+1);
fprintf("%f - %f\n",j1,j2)
oo = jag(:,(jag(1,:) >= j1) & ((jag(1,:)) <= j2));
nnn = size(oo,1);
[bc,be] = layered_bar_plot(oo, 2:nnn, 30, 316);

%% See? Random or not?
oo = jag(:,(jag(1,:) >= j1) & ((jag(1,:)) <= j2));
good = zeros(o1,o2);
for i = 1:size(oo,2)
    good(oo(2:mx+1,i), oo(mx+2:end,i)-50) = good(oo(2:mx+1,i), oo(mx+2:end,i)-50) +1;
end
figure(317);clf;
imagesc(good(1:fx,1:fx));
title("Detailed count of good hits");
colorbar();

%% z err as colormap
figure(318);clf;
zx = z(1:fx,1:fx);
x = zx - gt.gt_z(1:fx,1:fx);
x = abs(x);
imagesc(x);
title("z err as colormap (GT indication of inliers)");
colorbar();

%% gt.inl as colormap
figure(319);clf;
x = gt.inlmatrix(1:fx,1:fx);
imagesc(x);
title("gt.inl as colormap (GT inliers = 1)");
colorbar();

%% TODO; Just a save if you want to track this gt again
save("./" + "problematic_synth" + "_1_" + "gt.mat", 'gt');

%% Count the number of solutions in each try
o_low = -12;
o_hi  = mean(z(1:fx,1:fx),'all');
%o_hi= 12;

k = allsols(1,1:6*6);
k = allsols(1,:);
k = reshape(k,sx,[]);
% We're only interested in real solutions
% (though for numeric errors, we could perhaps allow some small
% imaginary part)
l = imag(k) == 0;
if 1
    % Filter out unreasonable o already?
    l2 = l;
    l2 = l2 & isfinite(k);
    l2 = l2 & (k <= o_hi);
    l2 = l2 & (k > o_low);
end
i = sum(l,1);
figure(300)
subplot(1,2,1);
histogram(i);
title("All real");
xlabel("Solutions from minsolver");
subplot(1,2,2);
i = sum(l2,1);
histogram(i);
title(sprintf("Filtered [%f,%f]",o_low,o_hi));
xlabel("Solutions from minsolver");

%%
for j = 0:sx
    x = l2;
    x(:,i ~= j) = false;
    x = x(:);
    fprintf("Giving %i useful solutions in %i cases total: %i \n", [j,size(k(:,i == j),2),nnz(x)]);
    if j == 0
        x = l2(:); % Show the unfiltered solutions in this case
    elseif ~any(x)
        continue; 
    end

    %% 
    x = allsols(1,x);
    nnn = mean(x(1,:));
    fprintf("avg selection: %.20e (gt %.20f diff %.20e )\n",nnn,gt.o(1),gt.o(1)-nnn)
    nnn = median(x(1,:));
    fprintf("med selection: %.20e (gt %.20f diff %.20e )\n",nnn,gt.o(1),gt.o(1)-nnn)
    nnn = min(x(1,:));
    fprintf("min selection: %.20e (gt %.20f diff %.20e )\n",nnn,gt.o(1),gt.o(1)-nnn)
    nnn = max(x(1,:));
    fprintf("max selection: %.20e (gt %.20f diff %.20e )\n",nnn,gt.o(1),gt.o(1)-nnn)

    figure(600+j)
    h = histogram(x,1e5);
    ylim([0,10])
    title("Peaks for " + num2str(j) + " solutions.")
    [a,b] = sort(h.Values);
    j1 = a(end-5:end);
    j2 = b(end-5:end);
    [h.BinEdges(j2)],[j1]
    j1 = a(end);
    j2 = b(end);
    [h.BinEdges(j2:j2+1)], [j1, size(jag)]

    x = sort(x);
    figure(700+j)
    plot(x);
    title("Flats for " + num2str(j) + " solutions.")
    %%
end

%% Can you use the solution groups directly?
    figure(599);
    x = k;
    x(~l2) = NaN;
    %x = sort(x); % sort each group
    x = x(:,1:1200);
    plot(x','.','MarkerSize',2);
    title("Solution groups directly layers wrong solutions away from correct?")

%% Is RANSAC doing a better job on this data?
iters = 10000;
tic;
[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type',type,'solver',oSolver,'display','iter','rank',oSolver.rank,'iters',iters);
t2 = toc;
fprintf("Done after %f s\n", t2)
check_offset_vector(bestsol,gt,'tol',Inf);

[a,b] = min(abs(stats2(4,:) - gt.o(1)));
fprintf('Closest %5d: inliers = %3d, loss = %e (o=%e, bias=%e, qual=%f)\n', stats2([1,2,3,4],b),a, stats2(2,b) - stats2(3,b));

% TODO: This strength value is perhaps better than current selection. 
% It doesn't take "one more inlier" if the loss is bad, but if loss is < 1,
% it is equivalent.
[a,b] = max(stats2(2,:) - stats2(3,:));
fprintf('Stg-est %5d: inliers = %3d, loss = %e (o=%e, bias=%e, qual=%f)\n', ...
    stats2([1,2,3,4],b), abs(stats2(4,b) - gt.o(1)), a);
