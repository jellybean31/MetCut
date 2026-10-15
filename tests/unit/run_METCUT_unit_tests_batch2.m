function R = run_METCUT_unit_tests_batch2()
%RUN_METCUT_UNIT_TESTS_BATCH2
% Second METCUT_v1 unit-test batch.
%
% Run from the METCUT_v1 root:
%
%   clear functions
%   rehash
%   addpath(fullfile(pwd,'tests','unit'));
%   R = run_METCUT_unit_tests_batch2;
%   rmpath(fullfile(pwd,'tests','unit'));

fprintf('\n============================================================\n');
fprintf('METCUT_v1 UNIT TEST SUITE - BATCH 2\n');
fprintf('============================================================\n');

rootDir = pwd;
functionDir = fullfile(rootDir,'functions');
unitDir = fullfile(rootDir,'tests','unit');

if ~isfolder(functionDir) || ~isfolder(unitDir)
    error('Run this suite from the METCUT_v1 root.');
end

originalPath = path;
guard = onCleanup(@() path(originalPath)); 
restoredefaultpath;
addpath(functionDir,'-begin');
addpath(unitDir,'-begin');
rehash;

tests = { ...
    'test_METCUT_effective_state_exactD_mix', ...
    'test_METCUT_subspace_principal_angles', ...
    'test_METCUT_synthetic_recovery_definitions'};

required = { ...
    'effective_state_exactD_mix', ...
    'subspace_principal_angles_deg', ...
    'synthetic_recovery_parameter_definitions', ...
    'synthetic_recovery_make_truth_sample'};

missing = {};
for i = 1:numel(required)
    if isempty(which(required{i}))
        missing{end+1} = required{i}; 
    end
end

if ~isempty(missing)
    fprintf('Missing required functions:\n');
    for i = 1:numel(missing)
        fprintf('  - %s\n',missing{i});
    end
    error('METCUT:Batch2:MissingFunctions', ...
        'Required functions are missing.');
end

results = struct('name',{},'pass',{},'seconds',{},'errorMessage',{});

for i = 1:numel(tests)
    t0 = tic;
    S = struct('name',tests{i},'pass',false,'seconds',NaN,'errorMessage','');

    fprintf('\n[%d/%d] %s\n',i,numel(tests),tests{i});
    fprintf('------------------------------------------------------------\n');

    try
        out = feval(tests{i});
        S.pass = isstruct(out) && isfield(out,'pass') && logical(out.pass);
        if ~S.pass
            error('%s did not return PASS.',tests{i});
        end
    catch ME
        S.errorMessage = ME.message;
    end

    S.seconds = toc(t0);
    results(end+1) = S; %#ok<AGROW>
end

nPass = nnz([results.pass]);
nFail = numel(results)-nPass;

fprintf('\n============================================================\n');
fprintf('METCUT_v1 UNIT TEST SUITE - BATCH 2 SUMMARY\n');
fprintf('============================================================\n');

for i = 1:numel(results)
    if results(i).pass
        status = 'PASS';
    else
        status = 'FAIL';
    end

    fprintf('%-52s %s  %.3f s\n', ...
        results(i).name,status,results(i).seconds);

    if ~results(i).pass
        fprintf('    %s\n',results(i).errorMessage);
    end
end

fprintf('------------------------------------------------------------\n');
fprintf('Passed: %d / %d\n',nPass,numel(results));
fprintf('Failed: %d\n',nFail);
fprintf('============================================================\n');

R = struct;
R.suite = 'METCUT_v1 unit tests batch 2';
R.results = results;
R.nPass = nPass;
R.nFail = nFail;
R.pass = nFail == 0;

assert(R.pass,'METCUT:Batch2:Failure', ...
    '%d unit test(s) failed.',nFail);
end
