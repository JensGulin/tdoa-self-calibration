function [sols,coeffs] = solver_cotoa_rank3_55(z)
% SOLVER_COTOA_RANK3_55 Solve TxOA offsets for a specific case.
%   sols = SOLVER_COTOA_RANK3_55(z) solves for the offset o, given the
%       matrix of TxOA measurements z. The measurements are
%       given by z_ij = || r_i - s_j || + o_j, where r_i and s_j are
%       (unknown) receiver and senders positions. 
% 
%       COTOA means that there is a constant o_j = o for all j (and i).
%       Rank 3 can handle up to 3D positions, thus also 2D.
%       There are multiple putative offsets given, approximations assuming
%       small gaussian errors and no outliers.
% Input: 
%     z - 5x5 matrix with time-difference-of-arrival
%         (COTOA) measurements of the true TOA d, affected by a 
%         common offset o, i.e. z(i,j) = d(i,j) + o + e(i,j)
%         where e may be a small additive normal noise.
%         
% Output:
%     sols - 4x1 vector with all putative offsets o
%         so that [z(i,j) - o] produce a matrix with 
%         time-of-arrival measurements. 

% From ICASSP 2019 "ROBUST SELF-CALIBRATION OF CONSTANT OFFSET
%   TIME-DIFFERENCE-OF-ARRIVAL"
% f(o) = det(C' * (z − o).^2 * C) == 0
% In this case we only want o.
% Expand (z-o).^2 as z.^2 -2z.*o + o.^2
% and transform each term A into C'*A*C to get (see code)
% f(o) = det(Cz2C + C2zC .* -o + 0) == 0
% (Keeping the sign with o instead of C2zC hints to disregard
% the sign now, and just negate the roots at the end.)
% Note that while o is a matrix of same size as z here, it's actually
% a constant distributed to all elements, thus C'(o.^2)C = 0
% (also, the element-wise multiplication can be separated from C2zC).

C = [-ones(1,4) ; eye(4)];
Cz2C = C'*(z.^2)*C;
C2zC = C'*(2*z)*C;

%% Different approach: for A-oB, take the eigenvalues of inv(B)*A,
% since the determinant is the product of eigenvalues.
% f(o) = det(Cz2C + C2zC .* -o + 0) == 0
% f(o) = det( C2zC * ( inv(C2zC)*Cz2C - I.*o) )
%      = det(C2zC) * prod( eig( inv(C2zC)*Cz2C ) - o ) == 0
% Any term of the product can be zero, exactly when eig == o.
% With det(C2zC) non-zero, C2zC is invertible and the eig provides the roots.
%    s = eig( inv(C2zC)*Cz2C );
% QZ algorithm (generalized Schur) skips the inverse, thus works for missing inverse.
% TODO: What cases need the generalized eigenvalues, are they more stable in any case?
% Note, det(C2zC) == 0 could mean that o can't be solved this way in any case.
% TODO: Could double check and discard where r < 3. 
% (two or more identical generalized eigenvalues, and
% corresponding eigenvectors span more than one independent direction.)
%  tol = 1e-9; r = rank(A + o(k)*B, tol)
    
s = eig( Cz2C,C2zC );
sols = ones(5,1)*s'; % Seems a bit faster than repmat.
return;

%%

% Setting up a template that indexes data to calculate the determinant
% (as coeffs [c4, c3, c2, c1, c0], with f(o) = c4*o^4 + c3*o^3 + c2*o^2 + c1*o + c0).
% Method is not clearly a standard Laplace expansion (for determinant), but
% all combinations are calculated through the intermediaries (minors), 
% prod1 and prod2, then further combined according to ids3. The sign and
% order is prepared by id1 and id3 respectively, so simple sums provides
% coeffs (in degree of o).

data = [Cz2C(:) ; C2zC(:)];
ids1 = [ ...
    1 , 1 , 17 , 17 , 3 , 3 , 19 , 19 , 5 , 21 , 5 , 21 , 6 , 22 , 6 , 22 , 1 , 1 , 17 , 17 , 2 , 2 , 18 , 18 ; ...
    6 , 22 , 6 , 22 , 8 , 24 , 8 , 24 , 3 , 3 , 19 , 19 , 4 , 4 , 20 , 20 , 8 , 24 , 8 , 24 , 7 , 23 , 7 , 23 ; ...
    5 , 21 , 5 , 21 , 7 , 23 , 7 , 23 , 1 , 1 , 17 , 17 , 2 , 2 , 18 , 18 , 5 , 21 , 5 , 21 , 6 , 22 , 6 , 22 ; ...
    2 , 2 , 18 , 18 , 4 , 4 , 20 , 20 , 7 , 23 , 7 , 23 , 8 , 24 , 8 , 24 , 4 , 4 , 20 , 20 , 3 , 3 , 19 , 19 ];

ids2 = [ ...
    9 , 9 , 25 , 25 , 11 , 11 , 27 , 27 , 9 , 9 , 25 , 25 , 10 , 10 , 26 , 26 , 9 , 9 , 25 , 25 , 10 , 10 , 26 , 26 ; ...
    14 , 30 , 14 , 30 , 16 , 32 , 16 , 32 , 15 , 31 , 15 , 31 , 16 , 32 , 16 , 32 , 16 , 32 , 16 , 32 , 15 , 31 , 15 , 31 ; ...
    13 , 29 , 13 , 29 , 15 , 31 , 15 , 31 , 13 , 29 , 13 , 29 , 14 , 30 , 14 , 30 , 13 , 29 , 13 , 29 , 14 , 30 , 14 , 30 ; ...
    10 , 10 , 26 , 26 , 12 , 12 , 28 , 28 , 11 , 11 , 27 , 27 , 12 , 12 , 28 , 28 , 12 , 12 , 28 , 28 , 11 , 11 , 27 , 27 ];

prod1 = data(ids1(1,:)).*data(ids1(2,:)) - ...
    data(ids1(3,:)).*data(ids1(4,:));
prod2 = data(ids2(1,:)).*data(ids2(2,:)) - ...
    data(ids2(3,:)).*data(ids2(4,:));

ids3 = [ ...
    4 ,  8 , 12 , 16 , 20 , 24 ,  2 ,  3 ,  4 ,  4 ,  6 ,  7 ,  8 ,  8 , 10 , 11 , 12 , 12 , 14 , 15 , 16 , 16 , 18 , 19 , 20 , 20 , 22 , 23 , 24 , 24 ,  1 ,  2 ,  2 ,  3 ,  3 ,  4 ,  5 ,  6 ,  6 ,  7 ,  7 ,  8 ,  9 , 10 , 10 , 11 , 11 , 12 , 13 , 14 , 14 , 15 , 15 , 16 , 17 , 18 , 18 , 19 , 19 , 20 , 21 , 22 , 22 , 23 , 23 , 24 ,  1 ,  1 ,  2 ,  3 ,  5 ,  5 ,  6 ,  7 ,  9 ,  9 , 10 , 11 , 13 , 13 , 14 , 15 , 17 , 17 , 18 , 19 , 21 , 21 , 22 , 23 ,  1 ,  5 ,  9 , 13 , 17 , 21  ; ...
    8 ,  4 , 16 , 12 , 24 , 20 ,  8 ,  8 ,  6 ,  7 ,  4 ,  4 ,  2 ,  3 , 16 , 16 , 14 , 15 , 12 , 12 , 10 , 11 , 24 , 24 , 22 , 23 , 20 , 20 , 18 , 19 ,  8 ,  6 ,  7 ,  6 ,  7 ,  5 ,  4 ,  2 ,  3 ,  2 ,  3 ,  1 , 16 , 14 , 15 , 14 , 15 , 13 , 12 , 10 , 11 , 10 , 11 ,  9 , 24 , 22 , 23 , 22 , 23 , 21 , 20 , 18 , 19 , 18 , 19 , 17 ,  6 ,  7 ,  5 ,  5 ,  2 ,  3 ,  1 ,  1 , 14 , 15 , 13 , 13 , 10 , 11 ,  9 ,  9 , 22 , 23 , 21 , 21 , 18 , 19 , 17 , 17 ,  5 ,  1 , 13 ,  9 , 21 , 17 ];
prod3 = prod1(ids3(1,:)).*prod2(ids3(2,:));

coeffs = [ sum(prod3(1:6)) sum(prod3(7:30)) sum(prod3(31:66)) sum(prod3(67:90)) sum(prod3(91:96)) ];

% Note: Negating roots at the end, as we solved for -o so far.
sols = -roots(coeffs);

sols = ones(5,1)*sols'; % Seems a bit faster than repmat.
return;

%% Debug
%{
ids1 make minors for col 1+2, ids2 for col 3+4. ids3 combine them.
The end result (when solving -o, not o) should be this:
= X11 X22 X33 X44 + X11 X23 X34 X42 + X11 X24 X32 X43
- X11 X24 X33 X42 - X11 X23 X32 X44 - X11 X22 X34 X43
- X12 X21 X33 X44 - X13 X21 X34 X42 - X14 X21 X32 X43
+ X14 X21 X33 X42 + X13 X21 X32 X44 + X12 X21 X34 X43
+ X12 X23 X31 X44 + X13 X24 X31 X42 + X14 X22 X31 X43
- X14 X23 X31 X42 - X13 X22 X31 X44 - X12 X24 X31 X43
- X12 X23 X34 X41 - X13 X24 X32 X41 - X14 X22 X33 X41
+ X14 X23 X32 X41 + X13 X22 X34 X41 + X12 X24 X33 X41
Hint: The sign is negative if there is an odd number of 
"inversions" to get the column numbers in order. 
E.g. X11 X24 X33 X42 => X11 (X33 X24) X42 =>
 X11 X33 (X42 X24) => X11 (X42 X33) X24
 are three inversions.

% The "no o" coeff comes entirely from Cz2C. (close enough?)
assert(det(Cz2C) == coeffs(end), "c0 error: " + num2str(det(Cz2C) - coeffs(end)));
% The "all o" coeff comes entirely from C2zC.
assert(det(C2zC) == coeffs(1),   "c4 error: " + num2str(det(C2zC) - coeffs(1)));

%% Indexes in data correspond to stacked cols of [Cz2C C2zC].
% When col > c, it's C2zC
[row, c] = size(Cz2C);
assert(row == c);
j = c + 1; % number of coeffs
[row1, col1] = ind2sub(size([Cz2C C2zC]), ids1);
col1 = mod(col1-1,c)+1; % See symmetry without o
[row2, col2] = ind2sub(size([Cz2C C2zC]), ids2);
col2 = mod(col2-1,c)+1; % See symmetry without o
% Note: ids3 index does not follow the same system

%% Arrange strings to help understand
[row, col] = ind2sub(size(Cz2C), 1:numel(Cz2C));
X = compose('(%d,%d)', row', col')';
datapos = ["A" + X , "-oB" + X ];
datapos = ["A" + X , "oB" + X ]; % Solve for -o and negate later.
datadeg = [zeros(size(datapos,2)/2,1) , ones(size(datapos,2)/2,1) ];

%% prod1
Y = @(a) [a(1,:)+a(2,:);a(3,:)+a(4,:)]';
signature = ["","--"]; % Highlight, with "--", minors that were negated
flipped = diff(col1([1,2],:))<0;
prod1pos = compose("%s{ %s%s - %s%s }", signature(flipped+1)', datapos(ids1)');
prod1deg = Y(datadeg(ids1));
assert(~any(diff(prod1deg'))); prod1deg(:,2) = [];

%% prod2
flipped = diff(col2([1,2],:))<0;
prod2pos = compose("%s{ %s%s - %s%s }", signature(flipped+1)', datapos(ids2)');
prod2deg = Y(datadeg(ids2));
assert(~any(diff(prod2deg'))); prod2deg(:,2) = [];

%% Some terms should be negative, unless already reversed in prod1
% not odd distance && not flipped
flipped = diff(col1([1,2],:))<0;
col = ~mod(diff(row1([1,2],:)),2)' .* ~flipped';
col = find(col);
ids1f = ids1;
% Negate, by swapping the terms around the minus sign.
ids1f(:, col) = ids1f([3 4 1 2],col);
[ids1(:,col); nan(size(col')); ids1f(:,col)]

%% prod3
prod3pos = compose("%s * %s", prod1pos(ids3(1,:)), prod2pos(ids3(2,:)));
prod3deg = prod1deg(ids3(1,:)) + prod2deg(ids3(2,:));

%% Check that sign is correct now
% Elements with even col diff should be negative.
selection = cell(1,j);
for i=1:j
    fprintf("** %d **\n",j-i);
    [row, col] = find(prod3deg' == j-i);
    assert(all(row == 1));
    assert(all(prod3deg(col) == j-i));
    selection(i) = {col};
    % Output: "idx", "Should be negated", "total prod3"
    compose("%2d %s ( %s )", col', signature(~mod(diff(row1([1,2],ids3(1,col))),2)+1)', prod3pos(col,:))
    assert(sum(prod3(col)) == coeffs(i))
end

%% Reorder ids3 to get the coeff ordering right
ids3f = ids3(:,[selection{:}])
assert(all([selection{:}] == 1:size(ids3,2)), "order is not right, update ids3 = ids3f");

%% Speed test (not correct output, but indicating timing)
% The final (current) version saves a few microseconds...
iter = 1e6;
tic;
for i=1:iter
CC2 = sparse([5 4 4 3 4 3 3 2 4 3 3 2 3 2 2 1 5 4 4 3 4 3 3 2 4 3 3 2 3 2 2 1 5 4 4 3 4 3 3 2 4 3 3 2 3 2 2 1 5 4 4 3 4 3 3 2 4 3 3 2 3 2 2 1 5 4 4 3 4 3 3 2 4 3 3 2 3 2 2 1 5 4 4 3 4 3 3 2 4 3 3 2 3 2 2 1],1:96,[ones(1,32) -ones(1,32) ones(1,32)],5,96);
coeffs = CC2*prod3;
end
toc;

%%
[row,col] = find(CC2');
tic;
for i=1:iter
CC3 = prod3(row);
coeffs = [ sum(CC3(1:6)) sum(CC3(6+(1:24))) sum(CC3(30+(1:36))) sum(CC3(66+(1:24))) sum(CC3(90+(1:6))) ]';
end
toc;

%%
tic;
for i=1:iter
coeffs = [ sum(prod3(1:6)) sum(prod3(6+(1:24))) sum(prod3(30+(1:36))) sum(prod3(66+(1:24))) sum(prod3(90+(1:6))) ]';
end
toc;

%% Final version
tic;
for i=1:iter
data = [Cz2C(:) ; C2zC(:)];
ids1 = [ ...
    1 , 1 , 17 , 17 , 3 , 3 , 19 , 19 , 5 , 21 , 5 , 21 , 6 , 22 , 6 , 22 , 1 , 1 , 17 , 17 , 2 , 2 , 18 , 18 ; ...
    6 , 22 , 6 , 22 , 8 , 24 , 8 , 24 , 3 , 3 , 19 , 19 , 4 , 4 , 20 , 20 , 8 , 24 , 8 , 24 , 7 , 23 , 7 , 23 ; ...
    5 , 21 , 5 , 21 , 7 , 23 , 7 , 23 , 1 , 1 , 17 , 17 , 2 , 2 , 18 , 18 , 5 , 21 , 5 , 21 , 6 , 22 , 6 , 22 ; ...
    2 , 2 , 18 , 18 , 4 , 4 , 20 , 20 , 7 , 23 , 7 , 23 , 8 , 24 , 8 , 24 , 4 , 4 , 20 , 20 , 3 , 3 , 19 , 19 ];

ids2 = [ ...
    9 , 9 , 25 , 25 , 11 , 11 , 27 , 27 , 9 , 9 , 25 , 25 , 10 , 10 , 26 , 26 , 9 , 9 , 25 , 25 , 10 , 10 , 26 , 26 ; ...
    14 , 30 , 14 , 30 , 16 , 32 , 16 , 32 , 15 , 31 , 15 , 31 , 16 , 32 , 16 , 32 , 16 , 32 , 16 , 32 , 15 , 31 , 15 , 31 ; ...
    13 , 29 , 13 , 29 , 15 , 31 , 15 , 31 , 13 , 29 , 13 , 29 , 14 , 30 , 14 , 30 , 13 , 29 , 13 , 29 , 14 , 30 , 14 , 30 ; ...
    10 , 10 , 26 , 26 , 12 , 12 , 28 , 28 , 11 , 11 , 27 , 27 , 12 , 12 , 28 , 28 , 12 , 12 , 28 , 28 , 11 , 11 , 27 , 27 ];

prod1 = data(ids1(1,:)).*data(ids1(2,:)) - ...
    data(ids1(3,:)).*data(ids1(4,:));
prod2 = data(ids2(1,:)).*data(ids2(2,:)) - ...
    data(ids2(3,:)).*data(ids2(4,:));

ids3 = [ ...
    4 ,  8 , 12 , 16 , 20 , 24 ,  2 ,  3 ,  4 ,  4 ,  6 ,  7 ,  8 ,  8 , 10 , 11 , 12 , 12 , 14 , 15 , 16 , 16 , 18 , 19 , 20 , 20 , 22 , 23 , 24 , 24 ,  1 ,  2 ,  2 ,  3 ,  3 ,  4 ,  5 ,  6 ,  6 ,  7 ,  7 ,  8 ,  9 , 10 , 10 , 11 , 11 , 12 , 13 , 14 , 14 , 15 , 15 , 16 , 17 , 18 , 18 , 19 , 19 , 20 , 21 , 22 , 22 , 23 , 23 , 24 ,  1 ,  1 ,  2 ,  3 ,  5 ,  5 ,  6 ,  7 ,  9 ,  9 , 10 , 11 , 13 , 13 , 14 , 15 , 17 , 17 , 18 , 19 , 21 , 21 , 22 , 23 ,  1 ,  5 ,  9 , 13 , 17 , 21  ; ...
    8 ,  4 , 16 , 12 , 24 , 20 ,  8 ,  8 ,  6 ,  7 ,  4 ,  4 ,  2 ,  3 , 16 , 16 , 14 , 15 , 12 , 12 , 10 , 11 , 24 , 24 , 22 , 23 , 20 , 20 , 18 , 19 ,  8 ,  6 ,  7 ,  6 ,  7 ,  5 ,  4 ,  2 ,  3 ,  2 ,  3 ,  1 , 16 , 14 , 15 , 14 , 15 , 13 , 12 , 10 , 11 , 10 , 11 ,  9 , 24 , 22 , 23 , 22 , 23 , 21 , 20 , 18 , 19 , 18 , 19 , 17 ,  6 ,  7 ,  5 ,  5 ,  2 ,  3 ,  1 ,  1 , 14 , 15 , 13 , 13 , 10 , 11 ,  9 ,  9 , 22 , 23 , 21 , 21 , 18 , 19 , 17 , 17 ,  5 ,  1 , 13 ,  9 , 21 , 17 ];
prod3 = prod1(ids3(1,:)).*prod2(ids3(2,:));

coeffs = [ sum(prod3(1:6)) sum(prod3(7:30)) sum(prod3(31:66)) sum(prod3(67:90)) sum(prod3(91:96)) ];

% Note: Negating roots at the end, as we solved for -o so far.
sols = -roots(coeffs);
end
toc;

%% Different approach, for A-oB, take the eigenvalues of inv(B)*A.
% f(o) = det(Cz2C + C2zC .* -o + 0) == 0
% f(o) = det( C2zC * ( inv(C2zC)*Cz2C - I*o) )
%      = det(C2zC) * prod( eig( inv(C2zC)*Cz2C ) - o ) == 0
% With det(C2zC) non-zero, C2zC is invertible and the eig provides the roots.
tic;
for i=1:iter
s = eig( inv(C2zC)*Cz2C );
end
toc;
%% With \ for inversion.
tic;
for i=1:iter
s = eig( C2zC \ Cz2C );
end
toc;
%% QZ algorithm (generalized Schur) skips the inverse and handles missing inverse.
tic;
for i=1:iter
s = eig( Cz2C,C2zC );
end
toc;

%}
