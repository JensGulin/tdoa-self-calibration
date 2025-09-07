function solver = get_offset_solver_func(func)
% GET_OFFSET_SOLVER_FUNC Setup the solver struct from the function.
%   solver = GET_OFFSET_SOLVER_FUNC(func) Fills the solver struct for m
%       by n matrix of TxOA measurements between m receivers and n senders,
%       of the specific dimensions rank, with the name and solv function.

    solver.solv = func;
    fname = split(func2str(solver.solv),'_');
    solver.m = str2double(fname{4}(end-1));
    solver.n = str2double(fname{4}(end));
    solver.rank = str2double(fname{3}(end));
    solver.name = upper(fname{2}) + ...
        sprintf(" Rank %i (%ir/%is)", ...
        solver.rank, solver.m, solver.n);
end
