function u=tdoa_calc_u_from_xyo(x,y,o)
% u=TDOA_CALC_U_FROM_XY(x,y,o)
% calculates distance measurements u from x and y
%  according to u_ij = || x_i - y_j || + o_ij
%  for all i to all j. When x and y are r and s, respectively,
%  and o is 1xn (TDOA) or constant or 1x1 (COTOA, TOA), 
%  u becomes z, but other uses are possible. 
% Output: 
%   u - mxn matrix
% Input: 
%   x - Mxm matrix
%   y - Nxn matrix
%   o - 1xn matrix when o is same within y_j and represents o_1j, or
%       mx1 matrix when o is same within x_i and represents o_i1 or
%       1x1 matrix when o is constant, representing o_11.
%  M and N could be any dimension really, but usually 2 or 3. 
% Note: When M and N differ, the effective dimension D = max(M,N),
%  and the smaller input is extended with zeros, as the best guess.

[x_dim,m] = size(x);
[y_dim,n] = size(y);

if x_dim > y_dim
    y((y_dim+1):x_dim,:)=zeros(x_dim-y_dim,n);
elseif y_dim > x_dim
    x((x_dim+1):y_dim,:)=zeros(y_dim-x_dim,n);
end

d = sqrt( sum( ( repmat(x,1,n) - repelem(y,1,m) ).^2 , 1 ) );
d = reshape(d,m,n);
%keyboard;
u = d+repmat(o,m,1);

%% Backstory: Picking the fastest implementation.
% Testing 100000 cycles [m,n] = [10,10]
% pdist2 : 2.797936 s
% repmat3: 0.473893 s
% pdist2 : 8.751291 s
% repmat : 0.435438 s
% repmat2: 0.600843 s
% bsxfun : 1.073961 s
% bsxfun2: 1.576046 s
% bsxfun3: 6.032687 s
% kron   : 1.395114 s

%{
cnt = 100000;
fprintf("Testing %d cycles [m,n] = [%d,%d]\n", cnt,m,n);
tic;
for i = 1:cnt
    g = pdist2(x',y');
end
t = toc;
fprintf("pdist2 : %f s\n",t);
assert(~any(g-g, 'all'));

tic;
for i = 1:cnt
d = ( ( vecnorm( repmat(x,1,n) - repelem(y,1,m) ) ) );
f = reshape(d,m,n);
end
t = toc;
fprintf("repmat3: %f s\n",t);
assert(~any(f-g>1e-15, 'all'));

tic;
for i = 1:cnt
%    e = pdist2(x',y');
%    e = pdist2(x',y', "fasteuclidean",CacheSize=10);
    e = pdist2(x',y', "fasteuclidean");
end
t = toc;
fprintf("pdist2 : %f s\n",t);
assert(~any(e-g>1e-15, 'all'));

tic;
for i = 1:cnt
d = sqrt( sum( ( repmat(x,1,n) - repelem(y,1,m) ).^2 , 1 ) );
e = reshape(d,m,n);
end
t = toc;
fprintf("repmat : %f s\n",t);
assert(~any(e-g, 'all'));

tic;
for i = 1:cnt
fun = @(a,b) sum((a - b).^2);
d = sqrt( ( fun( repmat(x,1,n), repelem(y,1,m) ) ) );
f = reshape(d,m,n);
end
t = toc;
fprintf("repmat2: %f s\n",t);
assert(~any(f-g, 'all'));

tic;
for i = 1:cnt
fun = @(a,b) sum((a - b).^2);
d = sqrt( ( bsxfun( fun, repmat(x,1,n), repelem(y,1,m) ) ) );
h = reshape(d,m,n);
end
t = toc;
fprintf("bsxfun : %f s\n",t);
assert(~any(h-g, 'all'));

tic;
for i = 1:cnt
% slow! j = table2array(combinations(1:5,1:5));
fun = @(a,b) sum((a - b).^2);
[j,d] = meshgrid(1:n,1:m);
d = sqrt( ( bsxfun( fun, x(:,d), y(:,j) ) ) );
j = reshape(d,m,n);
end
t = toc;
fprintf("bsxfun2: %f s\n",t);
assert(~any(j-g, 'all'));

tic;
for i = 1:cnt
fun = @(a,b) sum((x(:,a) - y(:,b)).^2)';
j = sqrt( ( bsxfun( fun, (1:m)', (1:n) ) ) );
end
t = toc;
fprintf("bsxfun3: %f s\n",t);
assert(~any(j-g, 'all'));

tic;
for i = 1:cnt
d = sqrt( sum( ( kron(ones(1,n),x) - kron(y,ones(1,m)) ).^2 , 1 ) );
d = reshape(d,m,n);
end
t = toc;
fprintf("kron   : %f s\n",t);
assert(~any(d-g, 'all'));

%}

end % function
