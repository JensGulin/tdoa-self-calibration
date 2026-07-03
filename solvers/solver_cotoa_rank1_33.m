function [sols,coeffs] = solver_cotoa_rank1_33(z)
% SOLVER_COTOA_RANK1_33 Solve TxOA offsets for a specific case.
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

C = [-ones(1,2) ; eye(2)];
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
sols = ones(4,1)*s'; % Seems a bit faster than repmat.
return;
