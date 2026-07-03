function sols = solver_tdoa_rank1_53(z)
% SOLVER_TDOA_RANK1_53 Solve TDOA offsets for the minimal case (5r,3s)
%   sols = SOLVER_TDOA_RANK2_74(z) solves for the TDOA offsets o, given an
%       n=7 by m=4 matrix of TDOA measurements z. The measurements are
%       given by z_ij = || r_i - s_j || + o_j, where r_i and s_j are
%       receiver and senders positions in 2D.

d2 = (z).^2;
u = [(d2(:,2:3)-repmat(d2(:,1),1,2)) (-2*z(:,2:end)) (2*z(:,1))]\ones(5,1);
sols = [u(end)/sum(u(1:2));u((1:2)+2)./u(1:2)];
