function [sols,coeffs] = solver_cotoa_rank2_44(z)
% SOLVER_COTOA_RANK2_44 Solve TxOA offsets for a specific case.
%   sols = SOLVER_COTOA_RANK2_44(z) solves for the offset o, given the
%       matrix of TxOA measurements z. The measurements are
%       given by z_ij = || r_i - s_j || + o_j, where r_i and s_j are
%       (unknown) receiver and senders positions. 
% 
%       COTOA means that there is a constant o_j = o for all j (and i).
%       Rank 2 can handle up to 2D positions.
%       There are multiple putative offsets given, approximations assuming
%       small gaussian errors and no outliers.
% Input: 
%     z - 4x4 matrix with time-difference-of-arrival
%         (COTOA) measurements of the true TOA d, affected by a 
%         common offset o, i.e. z(i,j) = d(i,j) + o + e(i,j)
%         where e may be a small additive normal noise.
%         
% Output:
%     sols - 3x1 vector with all putative offsets o
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

C = [-ones(1,3) ; eye(3)];
Cz2C = C'*(z.^2)*C;
C2zC = C'*(2*z)*C;

% Setting up a template that indexes data to calculate the determinant
% (as coeffs [c3, c2, c1, c0], with f(o) = c3*o^3 + c2*o^2 + c1*o + c0).
% Here col 3 is chosen for standard Laplace expansion (for determinant).
% Note: ids1 represents the minors, ids2 the coefficients
% of the corresponding weight. Then id3 combines the results. The sign and
% order is prepared by id1 and id3 respectively, so simple sums provides
% coeffs (in degree of o).

data = [Cz2C(:) ; C2zC(:)];

ids1 = [ ...
    1 , 1 , 10 , 10 , 2 , 2 , 11 , 11 , 4 , 13 , 4 , 13 ; ...
    5 , 14 , 5 , 14 , 6 , 15 , 6 , 15 , 3 , 3 , 12 , 12 ; ...
    4 , 13 , 4 , 13 , 5 , 14 , 5 , 14 , 1 , 1 , 10 , 10 ; ...
    2 , 2 , 11 , 11 , 3 , 3 , 12 , 12 , 6 , 15 , 6 , 15 ];
ids2 = [9 18 7 16 8 17];

prod1 = data(ids1(1,:)).*data(ids1(2,:)) - ...
    data(ids1(3,:)).*data(ids1(4,:));
prod2 = data(ids2);

ids3 = [ ...
    4 , 8 , 12 , 2 , 3 , 4 , 6 , 7 , 8 , 10 , 11 , 12 , 1 , 2 , 3 , 5 , 6 , 7 , 9 , 10 , 11 , 1 , 5 , 9 ; ...
    2 , 4 , 6 , 2 , 2 , 1 , 4 , 4 , 3 , 6 , 6 , 5 , 2 , 1 , 1 , 4 , 3 , 3 , 6 , 5 , 5 , 1 , 3 , 5 ];
prod3 = prod1(ids3(1,:)).*prod2(ids3(2,:));

coeffs = [ sum(prod3(1:3)) sum(prod3(4:12)) sum(prod3(13:21)) sum(prod3(22:24)) ];

% Note: Negating roots at the end, as we solved for -o so far.
sols = -roots(coeffs);

sols = ones(4,1)*sols'; % Seems a bit faster than repmat.
return;

%% Debug
%{

% The "no o" coeff comes entirely from Cz2C. (close enough?)
assert(det(Cz2C) == coeffs(end), "c0 error: " + num2str(det(Cz2C) - coeffs(end)));
% The "all o" coeff comes entirely from C2zC.
assert(det(C2zC) == coeffs(1),   "c3 error: " + num2str(det(C2zC) - coeffs(1)));

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
prod2pos = datapos(ids2)';
prod2deg = datadeg(ids2)';

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

%}
