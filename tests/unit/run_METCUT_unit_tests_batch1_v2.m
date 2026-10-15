function R = run_METCUT_unit_tests_batch1_v2()
fprintf('\n============================================================\n');
fprintf('METCUT_v1 UNIT TEST SUITE - BATCH 1 v2\n');
fprintf('============================================================\n');

rootDir=pwd;
functionDir=fullfile(rootDir,'functions');
unitDir=fullfile(rootDir,'tests','unit');
if ~isfolder(functionDir) || ~isfolder(unitDir)
    error('Run this suite from the METCUT_v1 root.');
end

originalPath=path;
guard=onCleanup(@()path(originalPath)); 
restoredefaultpath;
addpath(functionDir,'-begin');
addpath(unitDir,'-begin');
rehash;

tests={'test_METCUT_intracellular_water_exactD', ...
       'test_METCUT_reaction_incorporation_exactD', ...
       'test_METCUT_standing_wax_accumulation_exactD'};

required={'intracellular_water_dynamic_exactD', ...
          'wax_hydrogen_source_traits_postmalonate_null', ...
          'wax_hydrogen_source_pools_postmalonate_exactD', ...
          'wax_reaction_incorporation_traits_null', ...
          'wax_apply_reaction_incorporation_exactD', ...
          'alkane_biosynthesis_h_stoichiometry_reaction_explicit', ...
          'wax_combine_hydrogen_classes_exactD', ...
          'wax_accumulation_default_traits', ...
          'wax_accumulation_dynamic_exactD'};

missing={};
for i=1:numel(required)
    if isempty(which(required{i})), missing{end+1}=required{i}; end %#ok<AGROW>
end
if ~isempty(missing)
    fprintf('Missing required functions:\n');
    for i=1:numel(missing), fprintf('  - %s\n',missing{i}); end
    error('METCUT:Batch1v2:MissingFunctions','Required functions are missing.');
end

results=struct('name',{},'pass',{},'seconds',{},'errorMessage',{});
for i=1:numel(tests)
    t0=tic; S=struct('name',tests{i},'pass',false,'seconds',NaN,'errorMessage','');
    fprintf('\n[%d/%d] %s\n',i,numel(tests),tests{i});
    try
        out=feval(tests{i});
        S.pass=isstruct(out) && isfield(out,'pass') && logical(out.pass);
        if ~S.pass, error('%s did not return PASS.',tests{i}); end
    catch ME
        S.errorMessage=ME.message;
    end
    S.seconds=toc(t0);
    results(end+1)=S; %#ok<AGROW>
end

nPass=nnz([results.pass]); nFail=numel(results)-nPass;
fprintf('\n============================================================\n');
fprintf('BATCH 1 v2 SUMMARY\n');
fprintf('============================================================\n');
for i=1:numel(results)
    if results(i).pass, st='PASS'; else, st='FAIL'; end
    fprintf('%-52s %s  %.3f s\n',results(i).name,st,results(i).seconds);
    if ~results(i).pass, fprintf('    %s\n',results(i).errorMessage); end
end
fprintf('Passed: %d / %d\nFailed: %d\n',nPass,numel(results),nFail);

R=struct('results',results,'nPass',nPass,'nFail',nFail,'pass',nFail==0);
assert(R.pass,'METCUT:Batch1v2:Failure','%d unit test(s) failed.',nFail);
end
