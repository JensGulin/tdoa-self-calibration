function [zcalc,zok] = misstdoa_calc_z(sol)
% MISSTDOA_CALC_Z

m = length(sol.rows);
n = length(sol.cols);

switch lower(sol.type)
    case {'uvabo'}
        tmp = -2*(sol.u*sol.v) + repmat(sol.a,1,n) + repmat(sol.b,m,1);
        % This is (kind of) -2* U'*V + a_i + b_j
        ok = tmp >= 0;
        neg = tmp < 0;
        tmp = sqrt(relu(tmp)) + repmat(sol.o,m,1);
        % This is sqrt( -2*U'*V + a_i + b_j ) + o
        %ztmp = ztmp.*ok;
        
        %utmp = sqrt(abs(sol.w*sol.v+ repmat(sol.c_anchor,1,n)+repmat(sol.d2_anchor,m,1)))+repmat(sol.o,m,1);
        
        zcalc = nan(size(sol.z));
        zcalc(sol.rows,sol.cols) = tmp;
        
        zok = false(size(sol.z));
        zok(sol.rows,sol.cols) = ok;
    case {'rso'}
        tmp = toa_calc_d_from_xy(sol.r, sol.s);
        tmp = tmp + repmat(sol.o,m,1);
        
        zcalc = nan(size(sol.z));
        zcalc(sol.rows,sol.cols) = tmp;
        zok = sol.inlmatrix;
        zok = false(size(sol.z)); % TODO: all better than inl?
        zok(sol.rows,sol.cols) = true;
    case {'gt_rso'} % Ground truth from generate_synthetic_txoa()
        tmp = toa_calc_d_from_xy(sol.r, sol.s);
        zcalc = tmp + repmat(sol.o,m,1);
        zok = sol.inlmatrix;
    case {'z'} % Pre-solution illustration
    % TODO JAG Come up with a useful measure from z only?
        zcalc = repmat(mean(sol.z,'omitmissing'),m,1);
        zcalc = repmat(mean(sol.z,2,'omitmissing'),1,n);
        zok = sol.inlmatrix;
        %zok = sol.z;
    otherwise
        error("misstdoa_calc_z: Unknown sol type.");
end
