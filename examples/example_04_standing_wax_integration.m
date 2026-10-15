function R = example_04_standing_wax_integration()
%EXAMPLE_04_STANDING_WAX_INTEGRATION
% Run the integration example for standing plant-wax production,
% replacement, loss, and exact-D accumulation.


thisFile = mfilename('fullpath');
examplesDir = fileparts(thisFile);
root = fileparts(examplesDir);

fnDir = fullfile(root,'functions');
testDir = fullfile(root,'tests','integration');

originalPath = path;
cleanup = onCleanup(@() path(originalPath)); 
addpath(fnDir,'-begin');
addpath(testDir,'-begin');
rehash

fprintf('\n============================================================\n');
fprintf('METCUT EXAMPLE 04 - STANDING PLANT-WAX INTEGRATION\n');
fprintf('============================================================\n');

traits = wax_accumulation_default_traits();

fprintf('Default accumulation-trait fields:\n');
disp(fieldnames(traits))

fprintf('\nRunning validated standing-wax integration example...\n');
R = test_METCUT_standing_wax_integration;

fprintf('\nReturned result fields:\n');
disp(fieldnames(R))
end
