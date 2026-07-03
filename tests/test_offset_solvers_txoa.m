%TEST_OFFSET_SOLVERS_TXOA Test minimal offset solvers on synthetic TxOA data.
%
% This script checks that each listed minimal solver returns at least one
% offset vector matching the ground-truth offsets for noiseless, complete
% synthetic data. It is intentionally non-visual and assertion-based.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'components'));
addpath(fullfile(repoRoot, 'kalle_generalizations'));
addpath(fullfile(repoRoot, 'misc'));
addpath(fullfile(repoRoot, 'solvers'));
addpath(fullfile(repoRoot, 'systems'));
addpath(fullfile(repoRoot, 'tests'));
addpath(fullfile(repoRoot, 'util'));

cases = {
    'tdoa',  1, @solver_tdoa_rank1_53,  5, 3
    'tdoa',  2, @solver_tdoa_rank2_74,  7, 4
    'tdoa',  2, @solver_tdoa_rank2_56,  5, 6
    'tdoa',  3, @solver_tdoa_rank3_95,  9, 5
    'tdoa',  3, @solver_tdoa_rank3_76,  7, 6
    'tdoa',  3, @solver_tdoa_rank3_68,  6, 8
    'cotoa', 1, @solver_cotoa_rank1_33, 3, 3
    'cotoa', 2, @solver_cotoa_rank2_44, 4, 4
    'cotoa', 3, @solver_cotoa_rank3_55, 5, 5
    'toa',   1, @solver_toa_rank1_32,   3, 2
    'toa',   2, @solver_toa_rank2_43,   4, 3
    'toa',   3, @solver_toa_rank3_54,   5, 4
};

tol = 5e-4;
maxAttempts = 20;
seed0 = 1000;

fprintf('Running %d TxOA offset-solver tests...\n', size(cases, 1));

for caseIdx = 1:size(cases, 1)
    offsetType = cases{caseIdx, 1};
    rankValue = cases{caseIdx, 2};
    solverFunc = cases{caseIdx, 3};
    solverRows = cases{caseIdx, 4};
    solverCols = cases{caseIdx, 5};

    solver = get_offset_solver_func(solverFunc);
    assert(solver.rank == rankValue, ...
        'Solver %s reports rank %d, expected %d.', ...
        func2str(solverFunc), solver.rank, rankValue);
    assert(solver.m == solverRows && solver.n == solverCols, ...
        'Solver %s reports size %dx%d, expected %dx%d.', ...
        func2str(solverFunc), solver.m, solver.n, solverRows, solverCols);

    caseName = sprintf('%s rank %d %s', ...
        upper(offsetType), rankValue, func2str(solverFunc));

    [ok, usedSeed, err, lastMessage] = run_offset_solver_case( ...
        offsetType, rankValue, solverFunc, solverRows, solverCols, ...
        tol, seed0 + 100 * caseIdx, maxAttempts);

    assert(ok, ...
        'No valid offset solution found for %s after %d attempts. Last max error: %e. Last message: %s', ...
        caseName, maxAttempts, err, lastMessage);

    fprintf('  PASS %-42s seed=%d maxerr=%e\n', caseName, usedSeed, err);
end

fprintf('All TxOA offset-solver tests passed.\n');

function [ok, usedSeed, err, lastMessage] = run_offset_solver_case( ...
    offsetType, rankValue, solverFunc, solverRows, solverCols, ...
    tol, firstSeed, maxAttempts)

ok = false;
usedSeed = NaN;
err = Inf;
lastMessage = '';

for attemptIdx = 1:maxAttempts
    seed = firstSeed + attemptIdx - 1;

    [z, gt] = generate_synthetic_txoa( ...
        solverRows + 3, solverCols + 3, rankValue, ...
        offsetType, 0, 0, 0, [-2 6], seed);

    zMinimal = z(1:solverRows, 1:solverCols);
    gtOffset = gt.o(1:solverCols).';

    try
        sols = solverFunc(zMinimal);
    catch errObj
        lastMessage = errObj.message;
        continue
    end

    if size(sols, 1) ~= numel(gtOffset)
        lastMessage = sprintf( ...
            'Solver returned %d rows, expected %d.', ...
            size(sols, 1), numel(gtOffset));
        continue
    end

    solErr = max(abs(sols - gtOffset), [], 1);
    err = min(solErr);
    lastMessage = sprintf('best offset error was %e', err);
    if err < tol
        ok = true;
        usedSeed = seed;
        return
    end
end
end
