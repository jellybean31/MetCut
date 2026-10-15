function R = run_all_examples()
%RUN_ALL_EXAMPLES
% Run the five METCUT examples in sequence.


fprintf('\n============================================================\n');
fprintf('METCUT v1 - RUN ALL EXAMPLES\n');
fprintf('============================================================\n');

R = struct;

fprintf('\n[1/5] Reference forward run\n');
R.referenceForward = example_01_reference_forward_run();

fprintf('\n[2/5] Parameterized sensitivity\n');
R.parameterizedSensitivity = example_02_parameterized_sensitivity();

fprintf('\n[3/5] Precipitation tracer\n');
R.precipitationTracer = example_03_precipitation_tracer();

fprintf('\n[4/5] Standing plant-wax integration\n');
R.standingWax = example_04_standing_wax_integration();

fprintf('\n[5/5] Inverse recovery\n');
R.inverseRecovery = example_05_inverse_recovery();

fprintf('\n============================================================\n');
fprintf('ALL METCUT EXAMPLES COMPLETED\n');
fprintf('============================================================\n');
end
