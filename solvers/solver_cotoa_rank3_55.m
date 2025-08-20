function [sols,coeffs] = solver_cotoa_rank3_55(z)
% SOLVER_COTOA_RANK3_55 Solve TxOA offsets for a specific case.
%   sols = SOLVER_COTOA_RANK3_55(z) solves for the offset o, given the
%       matrix of TxOA measurements z. The measurements are
%       given by z_ij = || r_i - s_j || + o_j, where r_i and s_j are
%       (unknown) receiver and senders positions. 
% 
%       COTOA means that there is a constant o_j = o for all j (and i).
%       Rank 3 can handle up to 3D positions, but also 2D.
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
% f(o) = det(C' * (z − o).^2 * C) = 0
% In this case we only want o.
% Expand (z-o).^2 as z.^2 -2z.*o + o.^2
% and transform each term A into C'*A*C to get (see code)
% Cz2C + C2zC .* o + 0
% Note that while o is a matrix of same size as z here, it's actually
% a constant distributed to all elements, thus C'(o.^2)C = 0
% and the element-wise multiplication can be separated from C2zC.

C = [-ones(1,4) ; eye(4)];
Cz2C = C'*(z.^2)*C;
C2zC = C'*(2*z)*C; % Should be -, but get coeffs for -o and then -roots.
% Setting up a template that uses 
data = [Cz2C(:);C2zC(:)];
ids1 = [ ...
1 , 1 , 17 , 17 , 3 , 3 , 19 , 19 , 1 , 1 , 17 , 17 , 2 , 2 , 18 , 18 , 1 , 1 , 17 , 17 , 2 , 2 , 18 , 18 ; ...
6 , 22 , 6 , 22 , 8 , 24 , 8 , 24 , 7 , 23 , 7 , 23 , 8 , 24 , 8 , 24 , 8 , 24 , 8 , 24 , 7 , 23 , 7 , 23 ; ...
5 , 21 , 5 , 21 , 7 , 23 , 7 , 23 , 5 , 21 , 5 , 21 , 6 , 22 , 6 , 22 , 5 , 21 , 5 , 21 , 6 , 22 , 6 , 22 ; ...
2 , 2 , 18 , 18 , 4 , 4 , 20 , 20 , 3 , 3 , 19 , 19 , 4 , 4 , 20 , 20 , 4 , 4 , 20 , 20 , 3 , 3 , 19 , 19 ];
ids2 = [ ...
9 , 9 , 25 , 25 , 11 , 11 , 27 , 27 , 9 , 9 , 25 , 25 , 10 , 10 , 26 , 26 , 9 , 9 , 25 , 25 , 10 , 10 , 26 , 26 ; ...
14 , 30 , 14 , 30 , 16 , 32 , 16 , 32 , 15 , 31 , 15 , 31 , 16 , 32 , 16 , 32 , 16 , 32 , 16 , 32 , 15 , 31 , 15 , 31 ; ...
13 , 29 , 13 , 29 , 15 , 31 , 15 , 31 , 13 , 29 , 13 , 29 , 14 , 30 , 14 , 30 , 13 , 29 , 13 , 29 , 14 , 30 , 14 , 30 ; ...
10 , 10 , 26 , 26 , 12 , 12 , 28 , 28 , 11 , 11 , 27 , 27 , 12 , 12 , 28 , 28 , 12 , 12 , 28 , 28 , 11 , 11 , 27 , 27 ];
prod1 = data(ids1(1,:)).*data(ids1(2,:)) - ...
    data(ids1(3,:)).*data(ids1(4,:));
prod2 = data(ids2(1,:)).*data(ids2(2,:)) - ...
    data(ids2(3,:)).*data(ids2(4,:));
%
ids3 = [ ...
1 , 1 , 1 , 1 , 2 , 2 , 2 , 2 , 3 , 3 , 3 , 3 , 4 , 4 , 4 , 4 , 5 , 5 , 5 , 5 , 6 , 6 , 6 , 6 , 7 , 7 , 7 , 7 , 8 , 8 , 8 , 8 , 9 , 9 , 9 , 9 , 10 , 10 , 10 , 10 , 11 , 11 , 11 , 11 , 12 , 12 , 12 , 12 , 13 , 13 , 13 , 13 , 14 , 14 , 14 , 14 , 15 , 15 , 15 , 15 , 16 , 16 , 16 , 16 , 17 , 17 , 17 , 17 , 18 , 18 , 18 , 18 , 19 , 19 , 19 , 19 , 20 , 20 , 20 , 20 , 21 , 21 , 21 , 21 , 22 , 22 , 22 , 22 , 23 , 23 , 23 , 23 , 24 , 24 , 24 , 24 ; ...
5 , 6 , 7 , 8 , 5 , 6 , 7 , 8 , 5 , 6 , 7 , 8 , 5 , 6 , 7 , 8 , 1 , 2 , 3 , 4 , 1 , 2 , 3 , 4 , 1 , 2 , 3 , 4 , 1 , 2 , 3 , 4 , 13 , 14 , 15 , 16 , 13 , 14 , 15 , 16 , 13 , 14 , 15 , 16 , 13 , 14 , 15 , 16 , 9 , 10 , 11 , 12 , 9 , 10 , 11 , 12 , 9 , 10 , 11 , 12 , 9 , 10 , 11 , 12 , 21 , 22 , 23 , 24 , 21 , 22 , 23 , 24 , 21 , 22 , 23 , 24 , 21 , 22 , 23 , 24 , 17 , 18 , 19 , 20 , 17 , 18 , 19 , 20 , 17 , 18 , 19 , 20 , 17 , 18 , 19 , 20 ];
prod3 = prod1(ids3(1,:)).*prod2(ids3(2,:));

CC2 = sparse([5 4 4 3 4 3 3 2 4 3 3 2 3 2 2 1 5 4 4 3 4 3 3 2 4 3 3 2 3 2 2 1 5 4 4 3 4 3 3 2 4 3 3 2 3 2 2 1 5 4 4 3 4 3 3 2 4 3 3 2 3 2 2 1 5 4 4 3 4 3 3 2 4 3 3 2 3 2 2 1 5 4 4 3 4 3 3 2 4 3 3 2 3 2 2 1],1:96,[ones(1,32) -ones(1,32) ones(1,32)],5,96);
sols = -roots(CC2*prod3);
coeffs = CC2*prod3;

sols = ones(5,1)*sols';
