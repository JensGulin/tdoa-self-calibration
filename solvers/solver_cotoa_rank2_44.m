function [sols,coeffs] = tdoa_offset_44(a)
%sols = tdoa_offset_44(a)
% Input:
%     a - This is a 4x4 matrix a with time-difference-of-arrival
%         measurements
% Output:
%     sols - 3x1 vector with the three possible offsets o
%     so that a(i,j) - o produce a 4x4 matrix with
%     time-of-arrival measurements.

cl = [-ones(3,1) eye(3)];
cr = [-ones(3,1) eye(3)]';
c1 = cl*(a.^2)*cr;
cx = cl*(2*a)*cr;
data = [c1(:);cx(:)];
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
coeffs = [-sum(prod3(1:3)) sum(prod3(4:12)) -sum(prod3(13:20)) sum(prod3(21:24))];
sols = -roots(coeffs);

sols = ones(4,1)*sols';





