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
% f(o) = det(C' * (z − o).^2 * C) = 0
% In this case we only want o.
% Expand (z-o).^2 as z.^2 -2z.*o + o.^2
% and transform each term A into C'*A*C to get (see code)
% f(o) = det(Cz2C - C2zC .* o + 0) = 0
% Note that while o is a matrix of same size as z here, it's actually
% a constant distributed to all elements, thus C'(o.^2)C = 0
% and the element-wise multiplication can be separated from C2zC.

C = [-ones(1,3) ; eye(3)];
Cz2C = C'*(z.^2)*C;
C2zC = C'*(2*z)*C;

% Setting up a template that indexes data to calculate the determinant
% (as coeffs [c3, c2, c1, c0], with f(o) = c3*o^3 + c2*o^2 + c1*o + c0).
% Here col 3 is chosen for Laplace expansion (for determinant).
% Note: the - of C2zC is disregarded and brought in at the end.
% ids1 represents the minors, ids2 the coefficients
% of the corresponding weight. ids3 connects the output from
% ids1 and ids2, and reorders to group the coefficents nicely.

data = [Cz2C(:) ; C2zC(:)];

ids1 = [ ...
    1 , 1 , 10 , 10 , 2 , 2 , 11 , 11 , 3 , 3 , 12 , 12 ; ...
    5 , 14 , 5 , 14 , 6 , 15 , 6 , 15 , 4 , 13 , 4 , 13 ; ...
    4 , 13 , 4 , 13 , 5 , 14 , 5 , 14 , 6 , 15 , 6 , 15 ; ...
    2 , 2 , 11 , 11 , 3 , 3 , 12 , 12 , 1 , 1 , 10 , 10 ];
ids2 = [9 18 7 16 8 17];
prod1 = data(ids1(1,:)).*data(ids1(2,:)) - ...
    data(ids1(3,:)).*data(ids1(4,:));
prod2 = data(ids2);
ids3 = [ ...
    4 , 8 , 12 , 2 , 3 , 4 , 6 , 7 , 8 , 10 , 11 , 12 , 1 , 2 , 3 , 5 , 6 , 7 , 9 , 10 , 11 , 1 , 5 , 9 ; ...
    2 , 4 , 6 , 2 , 2 , 1 , 4 , 4 , 3 , 6 , 6 , 5 , 2 , 1 , 1 , 4 , 3 , 3 , 6 , 5 , 5 , 1 , 3 , 5 ];
prod3 = prod1(ids3(1,:)).*prod2(ids3(2,:));
coeffs = [-sum(prod3(1:3)) sum(prod3(4:12)) -sum(prod3(13:21)) sum(prod3(22:24))];
sols = roots(coeffs);

sols = ones(4,1)*sols';
return

%% Debug
%{
% The "no o" coeff comes entirely from Cz2C. (close enough?)
assert(det(Cz2C) == coeffs(end), "c0 error: " + num2str(det(Cz2C)- coeffs(end)));
% The "all o" coeff comes entirely from C2zC.
assert(-det(C2zC) == coeffs(1), "c3 error: " + num2str(-det(C2zC)- coeffs(1)));

% Indexes in data correspond to stacked cols of [Cz2C C2zC].
% When col > d, it's C2zC
[row1, col1] = ind2sub(size([Cz2C C2zC]), ids1);
[row2, col2] = ind2sub(size([Cz2C C2zC]), ids2);
% Note: ids3 index does not follow the same system

% Arrange strings to help understand
[row, col] = ind2sub(size(Cz2C), 1:numel(Cz2C));
X = arrayfun(@(a,b) sprintf('(%d,%d)', a, b), row, col, 'UniformOutput', false);
datapos = ["A" + X , "-oB" + X ];
prod1pos = datapos(ids1(1,:))' +  datapos(ids1(2,:))' +...
    repmat(' - ',12,1) +...
    datapos(ids1(3,:))' + datapos(ids1(4,:))';
prod2pos = datapos(ids2)';
% These terms should be negative, but already reversed in matching prod1
prod2pos(end-1:end) = ["--"; "--"] + prod2pos(end-1:end);

prod3pos = [prod1pos(ids3(1,:)) repmat('*',24,1) prod2pos(ids3(2,:))];
% odd o count means odd sign flips, compensate if A and B are without sign
prod3pos(1:3,:) % c3*o^3, odd o count, so needs minus
prod3pos(4:12,:) % c2*o^2, even o count
prod3pos(13:21,:) % c1*o^1, odd o count, so needs minus
prod3pos(22:24,:) % c0, even o count
%}
