function R = test_METCUT_synthetic_recovery_definitions()
%TEST_METCUT_SYNTHETIC_RECOVERY_DEFINITIONS


fprintf('\n============================================================\n');
fprintf('METCUT UNIT: SYNTHETIC-RECOVERY DEFINITIONS\n');
fprintf('============================================================\n');

D = synthetic_recovery_parameter_definitions();
details = struct;

assert(isstruct(D), ...
    'synthetic_recovery_parameter_definitions did not return a struct.');
assert(isfield(D,'transport') && isfield(D,'biochemistry') && isfield(D,'integrated'), ...
    'Synthetic-recovery registry is missing a required parameter group.');

details.nTransport = numel(D.transport);
details.nBiochemistry = numel(D.biochemistry);
details.nIntegrated = numel(D.integrated);

assert(details.nTransport == 2, ...
    'Synthetic-recovery transport registry should contain two parameters.');
assert(details.nBiochemistry == 7, ...
    'Synthetic-recovery biochemical registry should contain seven parameters.');
assert(details.nIntegrated == 9, ...
    'Integrated synthetic-recovery registry should contain nine parameters.');

% IDs must exist and be unique within the two component registries, without
% requiring the historical development labels themselves.
assert(isfield(D.transport,'id') && isfield(D.biochemistry,'id'), ...
    'Synthetic-recovery parameter definitions are missing id fields.');

transportIDs = string({D.transport.id});
biochemIDs = string({D.biochemistry.id});

assert(numel(unique(transportIDs)) == details.nTransport, ...
    'Transport parameter IDs are not unique.');
assert(numel(unique(biochemIDs)) == details.nBiochemistry, ...
    'Biochemical parameter IDs are not unique.');

% Deterministic truth-sample construction.
seed = 20260828;
T1 = synthetic_recovery_make_truth_sample(seed);
T2 = synthetic_recovery_make_truth_sample(seed);

details.truthWidth = width(T1);
details.truthHeight = height(T1);
details.sameSeedReproducible = isequaln(T1,T2);

assert(details.truthWidth == 10, ...
    'Truth-sample table should contain RunID plus nine parameters.');
assert(details.truthHeight >= 1, ...
    'Truth-sample constructor returned an empty table.');
assert(details.sameSeedReproducible, ...
    'Truth-sample construction is not reproducible for a fixed seed.');

R = struct;
R.name = 'test_METCUT_synthetic_recovery_definitions';
R.pass = true;
R.details = details;
R.truthSample = T1;

fprintf('Transport parameters:            %d\n',details.nTransport);
fprintf('Biochemical parameters:          %d\n',details.nBiochemistry);
fprintf('Integrated parameters:           %d\n',details.nIntegrated);
fprintf('Truth-sample width:              %d\n',details.truthWidth);
fprintf('Fixed-seed reproducible:         %d\n',details.sameSeedReproducible);
fprintf('PASS\n');
fprintf('============================================================\n');
end
