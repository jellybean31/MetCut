function R = example_05_inverse_recovery()
%EXAMPLE_05_INVERSE_RECOVERY
% Run the full-chain synthetic inverse-recovery example.
%
% The example compares final and time-resolved plant-wax observations.
% Time-resolved data improve L2 localization while
% I02 remains poorly constrained over the tested 2-20 d range.


thisFile = mfilename('fullpath');
examplesDir = fileparts(thisFile);
root = fileparts(examplesDir);

fnDir = fullfile(root,'functions');
testDir = fullfile(root,'tests','manuscript');

originalPath = path;
cleanup = onCleanup(@() path(originalPath)); 
addpath(fnDir,'-begin');
addpath(testDir,'-begin');
rehash

fprintf('\n============================================================\n');
fprintf('METCUT EXAMPLE 05 - FULL-CHAIN INVERSE RECOVERY\n');
fprintf('============================================================\n');
fprintf(['This example can take substantially longer than a single forward run\n' ...
         'because it evaluates a grid of candidate trajectories.\n\n']);

R = test_METCUT_manuscript_inverse_recovery;

fprintf('\nReturned result fields:\n');
disp(fieldnames(R))
end
