function R = run_METCUT_manuscript_tests()
%RUN_METCUT_MANUSCRIPT_TESTS


fprintf('\n============================================================\n');
fprintf('METCUT TEST SUITE\n');
fprintf('============================================================\n');

tests={ ...
    'test_METCUT_manuscript_precip_tracer', ...
    'test_METCUT_manuscript_source_memory', ...
    'test_METCUT_manuscript_wax_timing_targets', ...
    'test_METCUT_manuscript_inverse_recovery'};

n=numel(tests);
Name=strings(n,1);
Status=strings(n,1);
Elapsed_s=NaN(n,1);
Message=strings(n,1);
Output=cell(n,1);

for i=1:n
    Name(i)=string(tests{i});
    fprintf('\n[%d/%d] %s\n',i,n,Name(i));
    tt=tic;

    try
        Output{i}=feval(tests{i});
        Elapsed_s(i)=toc(tt);
        Status(i)="PASS";
        fprintf('PASS | %.3f s\n',Elapsed_s(i));
    catch ME
        Elapsed_s(i)=toc(tt);
        Status(i)="FAIL";

        if isempty(ME.identifier)
            Message(i)=string(ME.message);
        else
            Message(i)=string(ME.identifier)+": "+string(ME.message);
        end

        fprintf('FAIL | %.3f s\n%s\n',Elapsed_s(i),Message(i));
    end
end

Summary=table(Name,Status,Elapsed_s,Message);
disp(Summary);

R=struct;
R.version='2026-10-01_METCUT_manuscript_suite_final_v1';
R.summary=Summary;
R.output=Output;
R.nPass=sum(Status=="PASS");
R.nFail=sum(Status=="FAIL");
R.overallPass=R.nFail==0;

fprintf('Tests passed: %d / %d\n',R.nPass,n);
fprintf('Overall pass: %d\n',R.overallPass);
fprintf('============================================================\n');

assert(R.overallPass,'One or more METCUT tests failed.');
end
