function R = example_03_precipitation_tracer()
%EXAMPLE_03_PRECIPITATION_TRACER
% Run the precipitation-isotope tracer example used by the
% manuscript regression suite.

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
fprintf('METCUT EXAMPLE 03 - PRECIPITATION ISOTOPE TRACER\n');
fprintf('============================================================\n');
fprintf(['This example applies paired precipitation-isotope perturbations and\n' ...
         'tracks their propagation through the coupled isotope model.\n\n']);

R = test_METCUT_manuscript_precip_tracer;

fprintf('\nReturned result fields:\n');
disp(fieldnames(R))
end
