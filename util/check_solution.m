function check_solution(sol, gt, varargin)
%% Parse inputs.
p = inputParser;
valid_display = @(x) ismember(x, {'off', 'none', 'iter'});
addParameter(p, 'display', 'iter', valid_display);
addParameter(p, 'tol', 5e-4);
parse(p, varargin{:});
opts = p.Results;

print = opts.display == "iter";

%% 
plotsol = sol;
oo = gt.r;
oo = [oo; zeros(abs([3 0] - size(oo)))];
gt.r = oo;
oo = gt.s;
oo = [oo; zeros(abs([3 0] - size(oo)))];
gt.s = oo;
ooo = plotsol.r;
oo = zeros(3,size(ooo,2));
oo(1:size(ooo,1),plotsol.rows) = ooo;
r = oo;
ooo = plotsol.s;
oo = zeros(3,size(ooo,2));
oo(1:size(ooo,1),plotsol.cols) = ooo;
s = oo;
o = plotsol.o;
oo = zeros(1,size(ooo,2));
oo(1,plotsol.cols) = o;
o = oo;

%% Register result against ground truth (mic only)
if 0
[Q, t] = kabsch([r], [gt.r]);

r = Q * r + t;
s = Q * s + t;
end

%% Register result against ground truth.
[Q, t] = kabsch([r, s], [gt.r, gt.s]);

r = Q * r + t;
s = Q * s + t;


%% Plot results.
figure(2);
plot3(gt.r(1, :), gt.r(2, :), gt.r(3, :), 'ro');
hold on
plot3(gt.s(1, :), gt.s(2, :), gt.s(3, :), 'go');
plot3(r(1, :), r(2, :), r(3, :), 'r*');
plot3(s(1, :), s(2, :), s(3, :), 'g*');
axis equal
hold off

%% Print results.
er = vecnorm(r-gt.r);
es = vecnorm(s-gt.s);
eo = abs(o-gt.o);

fprintf('RMSE - r: %e, s: %e, o: %e\n', rms(er), rms(es), rms(eo));
fprintf('MAXE - r: %e, s: %e, o: %e\n', max(er), max(es), max(eo));
%% Check
if any([rms(er), rms(es), rms(eo)] > opts.tol)
    error("check_solution: Too large RMS!")
end
