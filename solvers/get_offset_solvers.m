function solvers = get_offset_solvers(varargin)
% GET_OFFSET_SOLVERS Get an array of solvers for finding TxOA offsets.
%   solvers = GET_OFFSET_SOLVERS() gets an array of solvers that take an m
%       by n matrix of TDOA measurements between m receivers and n senders,
%       and returns n offsets for the senders.
%   solvers = GET_OFFSET_SOLVERS(rank,m,n) filters the solvers by rank,
%       number of receivers or number of senders.
%   solvers = GET_OFFSET_SOLVERS(__,type) solvers for COTOA or TOA, instead
%       of the default TDOA. Can be an array, e.g. ["tdoa","cotoa"].
%
% Note: All variables are optional and can be named if out of order, e.g.
%   solvers = GET_OFFSET_SOLVERS('type','toa','m',5);

% Parse inputs.
p = inputParser;
valid_type = @(x) all(ismember(lower(x), {'tdoa', 'cotoa', 'toa'}));
addOptional(p, 'rank', []);
addOptional(p, 'm', []);
addOptional(p, 'n', []);
addOptional(p, 'type', "tdoa", valid_type);
parse(p, varargin{:});
opts = p.Results;
opts.type = lower(opts.type);

    solvers = [];
    
    if any(opts.type == "tdoa")

    % 2D 7/4 linear
    solvers(end+1).solv = @solver_tdoa_rank2_74;
    solvers(end).name = 'TDOA Rank 2 (7r/4s)';
    solvers(end).m = 7;
    solvers(end).n = 4;
    solvers(end).rank = 2;

    % 2D 5/6
    solvers(end+1).solv = @solver_tdoa_rank2_56;
    solvers(end).name = 'TDOA Rank 2 (5r/6s)';
    solvers(end).m = 5;
    solvers(end).n = 6;
    solvers(end).rank = 2;

    % 3D 9/5 linear
    solvers(end+1).solv = @solver_tdoa_rank3_95;
    solvers(end).name = 'TDOA Rank 3 (9r/5s)';
    solvers(end).m = 9;
    solvers(end).n = 5;
    solvers(end).rank = 3;

    % 3D 7/6
    solvers(end+1).solv = @solver_tdoa_rank3_76;
    solvers(end).name = 'TDOA Rank 3 (7r/6s)';
    solvers(end).m = 7;
    solvers(end).n = 6;
    solvers(end).rank = 3;

    % 3D 6/8
    solvers(end+1).solv = @solver_tdoa_rank3_68;
    solvers(end).name = 'TDOA Rank 3 (6r/8s)';
    solvers(end).m = 6;
    solvers(end).n = 8;
    solvers(end).rank = 3;

    end % tdoa
    
    if any(opts.type == "cotoa")
    
    solvers = [ solvers, ...
        get_offset_solver_func(@solver_cotoa_rank2_44), ... % 2D 4/4
        get_offset_solver_func(@solver_cotoa_rank3_55)  ... % 3D 5/5
        ];

    end % cotoa

    if any(opts.type == "toa")

    solvers = [ solvers, ...
        get_offset_solver_func(@solver_toa_rank2_43), ... % 2D 4/3
        get_offset_solver_func(@solver_toa_rank3_54)  ... % 3D 5/4
        ];

    end % toa

    if ~isempty(opts.rank)
        solvers([solvers.rank] ~= opts.rank) = [];
    end
    if ~isempty(opts.m)
        solvers([solvers.m] ~= opts.m) = [];
    end
    if ~isempty(opts.n)
        solvers([solvers.n] ~= opts.n) = [];
    end
end
