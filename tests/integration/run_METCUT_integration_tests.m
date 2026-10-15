function R = run_METCUT_integration_tests()
%RUN_METCUT_INTEGRATION_TESTS


fprintf('\n============================================================\n');
fprintf('METCUT INTEGRATION TEST SUITE\n');
fprintf('============================================================\n');

tests={ ...
    'test_METCUT_full_model_reference_regression', ...
    'test_METCUT_standing_wax_integration', ...
    'test_METCUT_hydraulic_solver_regression'};

n=numel(tests);
Name=strings(n,1);
Status=strings(n,1);
Elapsed_s=NaN(n,1);
Message=strings(n,1);
Output=cell(n,1);

for i=1:n
    Name(i)=string(tests{i});
    fprintf('\n[%d/%d] %s\n',i,n,Name(i));
    t=tic;
    try
        Output{i}=feval(tests{i});
        Elapsed_s(i)=toc(t);
        Status(i)="PASS";
        fprintf('PASS | %.3f s\n',Elapsed_s(i));
    catch ME
        Elapsed_s(i)=toc(t);
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
R.version='2026-09-30_METCUT_integration_suite_v2';
R.summary=Summary;
R.output=Output;
R.nPass=sum(Status=="PASS");
R.nFail=sum(Status=="FAIL");
R.overallPass=R.nFail==0;

fprintf('Integration tests passed: %d / %d\n',R.nPass,n);
fprintf('Overall integration pass: %d\n',R.overallPass);
fprintf('============================================================\n');

assert(R.overallPass,'One or more METCUT integration tests failed.');
end
