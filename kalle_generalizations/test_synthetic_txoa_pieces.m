%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test of individual solvers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (Find offset if no outliers)


%% CASE 1 - Generate synthetic data. TDOA RANK 3
m = 15;
n = 30;
dim = 3;
sigma = 1e-8; % Only a small gaussian noise
miss_ratio = 0.00;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'tdoa', sigma, miss_ratio, out_ratio);

% Test pieces

sols = solver_tdoa_rank3_95(z(1:9,1:5));
check_offset_vector(sols,gt.o(1:5)');
check_offset_vector(sols,gt.o(1:5)); % check accepts any transpose of o

sols = solver_tdoa_rank3_76(z(1:7,1:6));
check_offset_vector(sols,gt.o(1:6)');

sols = solver_tdoa_rank3_68(z(1:6,1:8));
check_offset_vector(sols,gt.o(1:8)');

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
check_offset_vector(sols,gt.o(1:4)');

sols = solver_tdoa_rank2_56(z(1:5,1:6));
check_offset_vector(sols,gt.o(1:6)');


%% CASE 2b - Generate synthetic data. TDOA RANK 1
m = 15;
n = 30;
dim = 1;
sigma = 1e-8;
miss_ratio = 0.00;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'tdoa', sigma, miss_ratio, out_ratio);

% Test pieces

sols = solver_tdoa_rank1_53a(z(1:5,1:3));
[sols gt.o(1:3)']
%check_offset_vector(sols,gt.o(1:3)');

figure(1)
hold off
plot(gt.r(1:5),zeros(1,5),'b*');
hold on
plot(gt.s(1:3),zeros(1,3),'bo');
for i = 1:3;
    text(gt.s(i),0.2,num2str(i),'FontSize',12);
end

% Det blir bara rätt för de ljud som är mellan mikrofonerna.

%%

if 0,
    addpath(genpath('/Users/kalle/Documents/projekt/autogen_v5'))

m = 5;
n = 3;
dim = 1;
sigma = 1e-8;
miss_ratio = 0.00;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'tdoa', sigma, miss_ratio, out_ratio);
    o = create_vars(3);

   tmp = compactionmatrix(5)*( (z-repmat(o',5,1)).^2 )*compactionmatrix(3)';

end



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
check_offset_vector(sols,gt.o(1:5)');

sols = solver_tdoa_rank3_76(z(1:7,1:6));
check_offset_vector(sols,gt.o(1:6)');

sols = solver_tdoa_rank3_68(z(1:6,1:8));
check_offset_vector(sols,gt.o(1:8)');

sols = solver_cotoa_rank3_55(z(1:5,1:5));
check_offset_vector(sols,gt.o(1:5)');


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
check_offset_vector(sols,gt.o(1:4)');

disp("This may fail: solver_tdoa_rank2_56, not all close enough")
sols = solver_tdoa_rank2_56(z(1:5,1:6));
check_offset_vector(sols,gt.o(1:6)','tol',1e-3);

sols = solver_cotoa_rank2_44(z(1:4,1:4));
check_offset_vector(sols,gt.o(1:4)');

% Can reuse rank3 solvers too!
sols = solver_cotoa_rank3_55(z(1:5,1:5));
check_offset_vector(sols,gt.o(1:5)');

%% CASE 4b - Generate synthetic data. COTOA RANK 1
m = 15;
n = 30;
dim = 1;
sigma = 1e-8;
miss_ratio = 0.00;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'cotoa', sigma, miss_ratio, out_ratio);

% Test pieces

sols = solver_tdoa_rank1_53(z(1:5,1:3));
check_offset_vector(sols,gt.o(1:3)');
% The 1D case is a bit special

% Can reuse rank3 solvers too!
sols = solver_cotoa_rank1_33(z(1:3,1:3));
check_offset_vector(sols,gt.o(1:3)');


%% CASE 5 - Generate synthetic data. TOA RANK 3
m = 15;
n = 30;
dim = 3;
sigma = 0;
miss_ratio = 0.00;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'toa', sigma, miss_ratio, out_ratio);

% Test pieces

sols = solver_tdoa_rank3_95(z(1:9,1:5));
check_offset_vector(sols,gt.o(1:5)');

sols = solver_tdoa_rank3_76(z(1:7,1:6));
check_offset_vector(sols,gt.o(1:6)');

sols = solver_tdoa_rank3_68(z(1:6,1:8));
check_offset_vector(sols,gt.o(1:8)');

sols = solver_cotoa_rank3_55(z(1:5,1:5));
check_offset_vector(sols,gt.o(1:5)');

sols = solver_toa_rank3_54(z(1:5,1:4));
check_offset_vector(sols,gt.o(1:4)');

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
check_offset_vector(sols,gt.o(1:4)');

sols = solver_tdoa_rank2_56(z(1:5,1:6));
check_offset_vector(sols,gt.o(1:6)');

sols = solver_cotoa_rank2_44(z(1:4,1:4));
check_offset_vector(sols,gt.o(1:4)');

sols = solver_toa_rank2_43(z(1:4,1:3));
check_offset_vector(sols,gt.o(1:3)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test of init_uvabo, extend, upgrade
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (Find full solution with missing and outliers)

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
asolver = get_offset_solver_func(@solver_tdoa_rank3_95);

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','tdoa','solver',asolver,'display','iter','rank',3,'iters',1000);
check_offset_vector(bestsol,gt);

[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')
solrso = upgrade_ransac(sol,'display', 'iter', 'threshold', 0.0001);

check_solution(solrso,gt);

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
asolver = get_offset_solver_func(@solver_tdoa_rank2_74);

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','tdoa','solver',asolver,'display','iter','rank',2,'iters',1000);
check_offset_vector(bestsol,gt);

[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')
solrso = upgrade_ransac(sol,'display', 'iter', 'threshold', 0.0001);

check_solution(solrso,gt);

%% CASE 2b - Generate synthetic data. TDOA RANK 1
m = 15;
n = 30;
dim = 1;
sigma = 1e-8;
miss_ratio = 0.02;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'tdoa', sigma, miss_ratio, out_ratio);

% Test pieces

% Default solver.
asolver = get_offset_solver_func(@solver_tdoa_rank1_53);

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','tdoa','solver',asolver,'display','iter','rank',1,'iters',1000);
check_offset_vector(bestsol,gt);

[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')
solrso = upgrade_ransac(sol,'display', 'iter', 'threshold', 0.0001);

check_solution(solrso,gt);



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
asolver = get_offset_solver_func(@solver_cotoa_rank3_55);

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','cotoa','solver',asolver,'display','iter','rank',3,'iters',1000);
check_offset_vector(bestsol,gt);

[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')
% TODO JAG Why refine only here?
[sol, res, jac] = refine_uvabo(sol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
solrso = upgrade_ransac(sol,'display', 'iter', 'threshold', 0.0001);

check_solution(solrso,gt);

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
asolver = get_offset_solver_func(@solver_cotoa_rank2_44);

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','cotoa','solver',asolver,'display','iter','rank',2,'iters',1000);
check_offset_vector(bestsol,gt);

[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')
solrso = upgrade_ransac(sol,'display', 'iter', 'threshold', 0.0001);

check_solution(solrso,gt);

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
asolver = get_offset_solver_func(@solver_toa_rank3_54);

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','toa','solver',asolver,'display','iter','rank',3,'iters',1000);
[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')
solrso = upgrade_ransac(sol,'display', 'iter', 'threshold', 0.0001);

check_solution(solrso,gt);


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
asolver = get_offset_solver_func(@solver_toa_rank2_43);

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','toa','solver',asolver,'display','iter','rank',2,'iters',1000);
[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')
solrso = upgrade_ransac(sol,'display', 'iter', 'threshold', 0.0001);

check_solution(solrso,gt);


%%










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Multi-dim options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test of init_uvabo, extend, upgrade
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (Find full solution with missing and outliers)

%% CASE 11 - Generate synthetic data. TDOA RANK 3-3
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
asolver = get_offset_solver_func(@solver_tdoa_rank3_95);

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','tdoa','solver',asolver,'display','iter','rank',3,'iters',1000);
[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')
solrso = upgrade_ransac(sol,'display', 'iter', 'threshold', 0.0001,'dims',dims);

check_solution(solrso,gt);

%% CASE 12 - Generate synthetic data. TDOA RANK 2-2
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
asolver = get_offset_solver_func(@solver_tdoa_rank2_74);

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','tdoa','solver',asolver,'display','iter','rank',2,'iters',1000);
[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')
solrso = upgrade_ransac(sol,'display', 'iter', 'threshold', 0.0001,'dims',dims);

check_solution(solrso,gt);

%% CASE 13 - Generate synthetic data. TDOA RANK 2-3
m = 15;
n = 30;
dims = [2 3];
dim = dims(1);
sigma = 1e-8;
miss_ratio = 0.02;
out_ratio = 0.00;

[z, gt] = generate_synthetic_txoa(m, n, dims, 'tdoa', sigma, miss_ratio, out_ratio);

% Test pieces

% Default solver.
asolver = get_offset_solver_func(@solver_tdoa_rank2_74);

[bestsol, max_inliers, best_err, stats1, stats2] = init_uvabo_ransac(z, 'offset_type','tdoa','solver',asolver,'display','iter','rank',dim,'iters',1000);
[sol, res, jac] = refine_uvabo(bestsol,'display','iter');
misstdoa_briefer_report(sol);
misstdoa_brief_visualization(sol);
sol = extend_uvabo_ransac(sol, 'display','iter')
solrso = upgrade_ransac(sol,'display', 'iter', 'threshold', 0.0001,'dims',dims);

check_solution(solrso,gt);

%% CASE 14 - Generate synthetic data. TDOA RANK 3-2
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
asolver = get_offset_solver_func(@solver_tdoa_rank2_74);
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

check_solution(solrso,gt);

%%






%% CASE 3 - Generate synthetic data. COTOA RANK 3
m = 15;
n = 30;
dim = [3 3];
sigma = 1e-8;
miss_ratio = 0.02;
out_ratio = 0.00;

% Generate data
[z, gt] = generate_synthetic_txoa(m, n, dim, 'cotoa', sigma, miss_ratio, out_ratio);

% Run system
[r, s, o, solrso] = txoa(z, 'display', 'iter', 'sigma', max(sigma, 1e-6), 'dims', [3 3], 'offset_type', 'cotoa');

% Check solution
check_solution(solrso,gt);






















%% TODO: Below is in test_syn_txoa.m also. Only there or only here?
% And it needs to be txoa generic.
% TODO: Check good thresholds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test full system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (Find full solution with missing and outliers)

%% CASE 1 - Generate synthetic data. TDOA RANK 3

%% Generate synthetic data.
m = 15;
n = 30;
dim = 3;
sigma = 1e-8;
sigma = 1e-4;
miss_ratio = 0.00;
out_ratio = 0.0;

[z, gt] = generate_synthetic_txoa(m, n, dim, 'tdoa', sigma, miss_ratio, out_ratio);

%% Run system.
[r, s, o, solrso] = tdoa(z, 'display', 'iter', 'sigma', max(sigma, 1e-6));
% [r, s, o, solrso] = tdoa_random(z, 'display', 'off', 'sigma', max(sigma, 1e-6), 'inits', 10);

check_solution(solrso,gt);
