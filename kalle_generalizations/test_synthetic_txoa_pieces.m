%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test of individual solvers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CASE 1 - Generate synthetic data. TDOA RANK 3
m = 15;
n = 30;
dim = 3;
sigma = 1e-8;
miss_ratio = 0.00;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'tdoa', sigma, miss_ratio, out_ratio);

% Test pieces

sols = solver_tdoa_rank3_95(z(1:9,1:5));
[sols gt.o(1:5)']

sols = solver_tdoa_rank3_76(z(1:7,1:6));
[sols gt.o(1:6)']

sols = solver_tdoa_rank3_68(z(1:6,1:8));
[sols gt.o(1:8)']

%% CASE 2 - Generate synthetic data. TDOA RANK 2
m = 15;
n = 30;
dim = 2;
sigma = 1e-8;
miss_ratio = 0.00;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'tdoa', sigma, miss_ratio, out_ratio);

% Test pieces

sols = solver_tdoa_rank2_74(z(1:7,1:4));
[sols gt.o(1:4)']

sols = solver_tdoa_rank2_56(z(1:5,1:6));
[sols gt.o(1:6)']

%% CASE 3 - Generate synthetic data. COTOA RANK 3
m = 15;
n = 30;
dim = 3;
sigma = 1e-8;
miss_ratio = 0.00;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'cotoa', sigma, miss_ratio, out_ratio);

% Test pieces

sols = solver_tdoa_rank3_95(z(1:9,1:5));
[sols gt.o(1:5)']

sols = solver_tdoa_rank3_76(z(1:7,1:6));
[sols gt.o(1:6)']

sols = solver_tdoa_rank3_68(z(1:6,1:8));
[sols gt.o(1:8)']

sols = solver_cotoa_rank3_55(z(1:5,1:5));
[sols gt.o(1:5)']


%% CASE 4 - Generate synthetic data. COTOA RANK 2
m = 15;
n = 30;
dim = 2;
sigma = 1e-8;
miss_ratio = 0.00;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'cotoa', sigma, miss_ratio, out_ratio);

% Test pieces

sols = solver_tdoa_rank2_74(z(1:7,1:4));
[sols gt.o(1:4)']

sols = solver_tdoa_rank2_56(z(1:5,1:6));
[sols gt.o(1:6)']

sols = solver_cotoa_rank2_44(z(1:4,1:4)); %Denna fungerar inte än
[sols gt.o(1:4)']

%% CASE 5 - Generate synthetic data. TOA RANK 3
m = 15;
n = 30;
dim = 3;
sigma = 1e-8;
miss_ratio = 0.00;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'toa', sigma, miss_ratio, out_ratio);

% Test pieces

sols = solver_tdoa_rank3_95(z(1:9,1:5));
[sols gt.o(1:5)']

sols = solver_tdoa_rank3_76(z(1:7,1:6));
[sols gt.o(1:6)']

sols = solver_tdoa_rank3_68(z(1:6,1:8));
[sols gt.o(1:8)']

sols = solver_cotoa_rank3_55(z(1:5,1:5));
[sols gt.o(1:5)']

sols = solver_toa_rank3_54(z(1:5,1:4));
[sols gt.o(1:4)']

%% CASE 6 - Generate synthetic data. TOA RANK 2
m = 15;
n = 30;
dim = 2;
sigma = 1e-8;
miss_ratio = 0.00;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'toa', sigma, miss_ratio, out_ratio);

% Test pieces

sols = solver_tdoa_rank2_74(z(1:7,1:4));
[sols gt.o(1:4)']

sols = solver_tdoa_rank2_56(z(1:5,1:6));
[sols gt.o(1:6)']

sols = solver_cotoa_rank2_44(z(1:4,1:4));
[sols gt.o(1:4)']

sols = solver_toa_rank2_43(z(1:4,1:3));
[sols gt.o(1:3)']


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test of init_uvabo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% CASE 1 - Generate synthetic data. TDOA RANK 3
m = 15;
n = 30;
dim = 3;
sigma = 1e-8;
miss_ratio = 0.01;
out_ratio = 0.01;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'tdoa', sigma, miss_ratio, out_ratio);

% Test pieces

% Default solver.
asolver.solv = @solver_tdoa_rank3_95;
asolver.name = func2str(asolver.solv);
asolver.m = 9;
asolver.n = 5;
asolver.rank = 3;

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','tdoa','solver',asolver,'display','iter','rank',3,'iters',1000);
[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')
solrso = upgrade_ransac(sol,'display', 'iter', 'threshold', 0.0001);


%% CASE 2 - Generate synthetic data. TDOA RANK 2
m = 15;
n = 30;
dim = 2;
sigma = 1e-8;
miss_ratio = 0.02;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'tdoa', sigma, miss_ratio, out_ratio);

% Test pieces

% Default solver.
asolver.solv = @solver_tdoa_rank2_74;
asolver.name = func2str(asolver.solv);
asolver.m = 7;
asolver.n = 4;
asolver.rank = 2;

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','tdoa','solver',asolver,'display','iter','rank',2,'iters',1000);
[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')
solrso = upgrade_ransac(sol,'display', 'iter', 'threshold', 0.0001);



%% CASE 3 - Generate synthetic data. COTOA RANK 3
m = 15;
n = 30;
dim = 3;
sigma = 1e-8;
miss_ratio = 0.02;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'cotoa', sigma, miss_ratio, out_ratio);

% Test pieces

% Default solver.
asolver.solv = @solver_cotoa_rank3_55;
asolver.name = func2str(asolver.solv);
asolver.m = 5;
asolver.n = 5;
asolver.rank = 3;

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','cotoa','solver',asolver,'display','iter','rank',3,'iters',1000);
[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
% TODO JAG Why only here?
sol = extend_uvabo_ransac(sol, 'display','iter')
[sol, res, jac] = refine_uvabo(sol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
solrso = upgrade_ransac(sol,'display', 'iter', 'threshold', 0.0001);


%% CASE 4 - Generate synthetic data. COTOA RANK 2
m = 15;
n = 30;
dim = 2;
sigma = 1e-8;
miss_ratio = 0.02;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'cotoa', sigma, miss_ratio, out_ratio);

% Test pieces

% Default solver.
asolver.solv = @solver_cotoa_rank2_44;
asolver.name = func2str(asolver.solv);
asolver.m = 4;
asolver.n = 4;
asolver.rank = 2;

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','cotoa','solver',asolver,'display','iter','rank',2,'iters',1000);
[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')



%% CASE 5 - Generate synthetic data. TOA RANK 3
m = 15;
n = 30;
dim = 3;
sigma = 1e-8;
miss_ratio = 0.02;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'toa', sigma, miss_ratio, out_ratio);
% Test pieces

% Default solver.
asolver.solv = @solver_toa_rank3_54;
asolver.name = func2str(asolver.solv);
asolver.m = 5;
asolver.n = 4;
asolver.rank = 3;

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','toa','solver',asolver,'display','iter','rank',3,'iters',1000);
[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')


%% CASE 6 - Generate synthetic data. TOA RANK 2
m = 15;
n = 30;
dim = 2;
sigma = 1e-8;
miss_ratio = 0.02;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'toa', sigma, miss_ratio, out_ratio);

% Test pieces

% Default solver.
asolver.solv = @solver_toa_rank2_43;
asolver.name = func2str(asolver.solv);
asolver.m = 4;
asolver.n = 3;
asolver.rank = 2;

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','toa','solver',asolver,'display','iter','rank',2,'iters',1000);
[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')


%%












%% CASE 11 - Generate synthetic data. TDOA RANK 3
m = 15;
n = 30;
dim = 3;
dims = [3 3];
sigma = 1e-8;
miss_ratio = 0.01;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'tdoa', sigma, miss_ratio, out_ratio);

% Test pieces

% Default solver.
asolver.solv = @solver_tdoa_rank3_95;
asolver.name = func2str(asolver.solv);
asolver.m = 9;
asolver.n = 5;
asolver.rank = 3;


[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','tdoa','solver',asolver,'display','iter','rank',3,'iters',1000);
[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')
solrso = upgrade_ransac(sol,'display', 'iter', 'threshold', 0.0001,'dims',dims);


%% CASE 12 - Generate synthetic data. TDOA RANK 2
m = 15;
n = 30;
dim = 2;
dims = [2 2];
sigma = 1e-8;
miss_ratio = 0.02;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'tdoa', sigma, miss_ratio, out_ratio);

% Test pieces

% Default solver.
asolver.solv = @solver_tdoa_rank2_74;
asolver.name = func2str(asolver.solv);
asolver.m = 7;
asolver.n = 4;
asolver.rank = 2;

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','tdoa','solver',asolver,'display','iter','rank',2,'iters',1000);
[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')
solrso = upgrade_ransac(sol,'display', 'iter', 'threshold', 0.0001,'dims',dims);


%% CASE 13 - Generate synthetic data. TDOA RANK 2
m = 15;
n = 30;
dim = [2 3];
sigma = 1e-8;
miss_ratio = 0.02;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'tdoa', sigma, miss_ratio, out_ratio);

% Test pieces

% Default solver.
asolver.solv = @solver_tdoa_rank2_74;
asolver.name = func2str(asolver.solv);
asolver.m = 7;
asolver.n = 4;
asolver.rank = 2;

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','tdoa','solver',asolver,'display','iter','rank',2,'iters',1000);
[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')
solrso = upgrade_ransac(sol,'display', 'iter', 'threshold', 0.0001,'dims',dim);


%% CASE 14 - Generate synthetic data. TDOA RANK 3 2
m = 15;
n = 30;
dim = [3 2];
dims = [3 2];
sigma = 1e-8;
miss_ratio = 0.02;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'tdoa', sigma, miss_ratio, out_ratio);

% Test pieces

% Default solver.
asolver.solv = @solver_tdoa_rank2_74;
asolver.name = func2str(asolver.solv);
asolver.m = 7;
asolver.n = 4;
asolver.rank = 2;

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','tdoa','solver',asolver,'display','iter','rank',2,'iters',1000);
[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')
solrso = upgrade_ransac(sol,'display', 'iter', 'threshold', 0.0001,'dims',dim);



%%































%% Run system.
[r, s, o, sol] = tdoa(z, 'display', 'iter', 'sigma', max(sigma, 1e-6));
% [r, s, o, sol] = tdoa_random(z, 'display', 'off', 'sigma', max(sigma, 1e-6), 'inits', 10);

%% Register result against ground truth.
[Q, t] = kabsch([r, s], [gt.r, gt.s]);

r = Q * r + t;
s = Q * s + t;

%% Plot results.
figure(2);
plot3(gt.r(1, :), gt.r(2, :), gt.r(3, :), 'r*');
hold on
plot3(gt.s(1, :), gt.s(2, :), gt.s(3, :), 'g*');
plot3(r(1, :), r(2, :), r(3, :), 'ro');
plot3(s(1, :), s(2, :), s(3, :), 'go');
axis equal
hold off

%% Print results.
er = vecnorm(r-gt.r);
es = vecnorm(s-gt.s);
eo = abs(o-gt.o);

fprintf('RMSE - r: %e, s: %e, o: %e\n', rms(er), rms(es), rms(eo));
fprintf('MAXE - r: %e, s: %e, o: %e\n', max(er), max(es), max(eo));
