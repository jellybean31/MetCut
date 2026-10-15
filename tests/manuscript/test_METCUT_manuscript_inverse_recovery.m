function R = test_METCUT_manuscript_inverse_recovery()
%TEST_METCUT_MANUSCRIPT_INVERSE_RECOVERY
%
% Regression for the full-chain wax-only synthetic
% inverse-recovery experiment.
%
% Fixture:
%   reconstructed from the package reference model
%   plus the synthetic observation:
%       sigmaWax_permil = 2.0
%       homologs = [27 29 31 33]


thisFile=mfilename('fullpath');
manuscriptDir=fileparts(thisFile);
testsDir=fileparts(manuscriptDir);
rootDir=fileparts(testsDir);
fnDir=fullfile(rootDir,'functions');
refDir=fullfile(rootDir,'data','reference');

originalPath=path;
cleanup=onCleanup(@() path(originalPath)); 
addpath(fnDir,'-begin');
rehash

hydPath=fullfile(refDir, ...
    'hydraulic_checkpoint_76cell_8day_hourly_corrected_mass_conservative_storagefixed.mat');
assert(isfile(hydPath),'Package-local hydraulic reference missing.');

%% 1. Transparent release fixture
[Mref,~]=run_full_model_from_hydraulic_checkpoint(hydPath);

Rpilot=struct;
Rpilot.truth=struct;
Rpilot.truth.model=Mref;
Rpilot.truth.config=Mref.config;

if ~isfield(Rpilot.truth.config,'biochemistry') || ...
        isempty(Rpilot.truth.config.biochemistry)
    Rpilot.truth.config.biochemistry=struct;
end

Rpilot.truth.config.biochemistry.sigmaWax_permil=2.0;
Rpilot.truth.config.biochemistry.homologs=[27 29 31 33];

%% 2. Accepted full-chain recovery design
options=struct;
options.I01Fixed=2.03690139;
options.truthI02=5.38306961;
options.nI02=11;
options.nL1=21;
options.nL2=21;
options.seed=20260828;

Q=full_chain_wax_only_grid_recovery(Rpilot,options);

%% 3. Structural release contract
assert(strcmp(Q.version,'2026-08-28_full_chain_wax_only_inverse_v1'), ...
    'Full-chain recovery version regression.');

assert(numel(Q.I02Grid)==11 && numel(Q.L1Grid)==21 && numel(Q.L2Grid)==21, ...
    'Full-chain recovery grid-size regression.');

assert(all(Q.passI02), ...
    'One or more upstream I02 candidates failed.');

assert(any(abs(double(Q.I02Grid(:))-options.truthI02)<=1e-12), ...
    'Exact truth I02 is missing from the recovery grid.');

assert(isequal(size(Q.final.J),[11 21 21]), ...
    'Final-only objective-grid shape regression.');
assert(isequal(size(Q.timeResolved.J),[11 21 21]), ...
    'Time-resolved objective-grid shape regression.');

nFiniteFinal=sum(isfinite(Q.final.J(:)));
nFiniteTime=sum(isfinite(Q.timeResolved.J(:)));

assert(nFiniteFinal==1925, ...
    'Final-only finite-objective-node count regression.');
assert(nFiniteTime==1925, ...
    'Time-resolved finite-objective-node count regression.');

assert(isequal(double(Q.truthFinal.homologs(:).'),[27 29 31 33]), ...
    'Final truth homolog contract regression.');
assert(isequal(double(Q.truthTimeResolved.homologs(:).'),[27 29 31 33]), ...
    'Time-resolved truth homolog contract regression.');
assert(isequaln(double(Q.truthFinal.sigma_permil),2.0), ...
    'Final truth observation-scale regression.');
assert(isequaln(double(Q.truthTimeResolved.sigma_permil),2.0), ...
    'Time-resolved truth observation-scale regression.');

assert(numel(Q.truthFinal.sampleIndex)==1, ...
    'Final-only observation count regression.');
assert(numel(Q.truthTimeResolved.sampleIndex)==17, ...
    'Time-resolved sample-count regression.');

%% 4. Accepted numerical fingerprints
truthEffective=effective5_truth_vector(Rpilot.truth);
truthLatent=effective5_to_latent2(truthEffective,Q.basis);

expectedTruth=[ ...
    5.38306961; ...
   -0.185074868056075; ...
   -0.0235422032918458];

assert(max(abs([options.truthI02;truthLatent(:)]-expectedTruth))<=1e-12, ...
    'Full-chain truth-coordinate regression.');

expectedFinalBest=[ ...
    20; ...
   -0.0976382208319103; ...
   -0.559192618314243];

expectedTimeBest=[ ...
    2; ...
   -0.195276441663821; ...
   -0.0931987697190405];

assert(max(abs(double(Q.final.best(:))-expectedFinalBest))<=1e-12, ...
    'Final-only best-fit coordinate regression.');
assert(max(abs(double(Q.timeResolved.best(:))-expectedTimeBest))<=1e-12, ...
    'Time-resolved best-fit coordinate regression.');

assert(abs(double(Q.final.Jmin)-7.57454832224391)<=1e-9, ...
    'Final-only Jmin regression.');
assert(abs(double(Q.timeResolved.Jmin)-689.007617455767)<=1e-8, ...
    'Time-resolved Jmin regression.');

assert(abs(double(Q.final.truthGridDeltaJ)-6.70680254154975)<=1e-9, ...
    'Final-only truth-grid deltaJ regression.');
assert(abs(double(Q.timeResolved.truthGridDeltaJ)-1.54656164200003)<=1e-9, ...
    'Time-resolved truth-grid deltaJ regression.');

%% 5. Manuscript-facing recovery behavior
assert(~logical(Q.final.truthInside68), ...
    'Final-only truth unexpectedly moved inside the joint 68%% region.');
assert(logical(Q.final.truthInside95), ...
    'Final-only truth is no longer inside the joint 95%% region.');

assert(logical(Q.timeResolved.truthInside68), ...
    'Time-resolved truth is no longer inside the joint 68%% region.');
assert(logical(Q.timeResolved.truthInside95), ...
    'Time-resolved truth is no longer inside the joint 95%% region.');

% I02 remains unconstrained across the complete tested domain.
I02=double(Q.I02Grid(:));
assertVectorEqualLocal(Q.final.I02Profile95,I02,1e-12, ...
    'Final-only I02 95%% profile regression.');
assertVectorEqualLocal(Q.timeResolved.I02Profile95,I02,1e-12, ...
    'Time-resolved I02 95%% profile regression.');

assert(abs(min(I02)-2)<=1e-12 && abs(max(I02)-20)<=1e-12, ...
    'I02 tested-domain regression.');

% L2 profile narrowing is the principal information-gain fingerprint.
finalL2_95=double(Q.final.L2Profile95(:));
timeL2_95=double(Q.timeResolved.L2Profile95(:));

expectedFinalL2_95=[ ...
   -0.559192618314243; ...
   -0.465993848595202; ...
   -0.372795078876162; ...
   -0.279596309157121; ...
   -0.186397539438081; ...
   -0.0931987697190405; ...
    0; ...
    0.0931987697190405];

expectedTimeL2_95=[ ...
   -0.279596309157121; ...
   -0.186397539438081; ...
   -0.0931987697190405; ...
    0];

assertVectorEqualLocal(finalL2_95,expectedFinalL2_95,1e-12, ...
    'Final-only L2 95%% profile regression.');
assertVectorEqualLocal(timeL2_95,expectedTimeL2_95,1e-12, ...
    'Time-resolved L2 95%% profile regression.');

finalL2Span=max(finalL2_95)-min(finalL2_95);
timeL2Span=max(timeL2_95)-min(timeL2_95);
spanRatio=timeL2Span/finalL2Span;

assert(abs(finalL2Span-0.652391388033283)<=1e-12, ...
    'Final-only L2 95%% span regression.');
assert(abs(timeL2Span-0.279596309157121)<=1e-12, ...
    'Time-resolved L2 95%% span regression.');
assert(abs(spanRatio-0.428571428571429)<=1e-12, ...
    'Time/final L2 95%% span-ratio regression.');
assert(timeL2Span<finalL2Span, ...
    'Time-resolved observations no longer narrow the L2 95%% profile.');

assert(numel(Q.final.L1Profile95)==1, ...
    'Final-only L1 95%% profile-width regression.');
assert(numel(Q.timeResolved.L1Profile95)==1, ...
    'Time-resolved L1 95%% profile-width regression.');

fprintf('\n');
fprintf('============================================================\n');
fprintf('METCUT FULL-CHAIN INVERSE-RECOVERY REGRESSION\n');
fprintf('============================================================\n');
fprintf('Upstream I02 candidates passed:     %d / %d\n', ...
    sum(Q.passI02),numel(Q.passI02));
fprintf('Finite objective nodes final/time:  %d / %d\n', ...
    nFiniteFinal,nFiniteTime);
fprintf('Truth [I02 L1 L2]:                  %s\n', ...
    mat2str(expectedTruth.',15));
fprintf('Final best [I02 L1 L2]:             %s\n', ...
    mat2str(double(Q.final.best(:).'),15));
fprintf('Time best [I02 L1 L2]:              %s\n', ...
    mat2str(double(Q.timeResolved.best(:).'),15));
fprintf('Truth inside 68%% final/time:        %d / %d\n', ...
    Q.final.truthInside68,Q.timeResolved.truthInside68);
fprintf('Truth inside 95%% final/time:        %d / %d\n', ...
    Q.final.truthInside95,Q.timeResolved.truthInside95);
fprintf('I02 95%% span final/time:            %.12g / %.12g d\n', ...
    max(Q.final.I02Profile95)-min(Q.final.I02Profile95), ...
    max(Q.timeResolved.I02Profile95)-min(Q.timeResolved.I02Profile95));
fprintf('L2 95%% span final/time:             %.15g / %.15g\n', ...
    finalL2Span,timeL2Span);
fprintf('L2 time/final span ratio:           %.15g\n',spanRatio);
fprintf('PASS\n');
fprintf('============================================================\n');

R=struct;
R.pass=true;
R.recovery=Q;
R.truthCoordinates=expectedTruth;
R.finalBest=double(Q.final.best(:));
R.timeResolvedBest=double(Q.timeResolved.best(:));
R.finalTruthGridDeltaJ=double(Q.final.truthGridDeltaJ);
R.timeResolvedTruthGridDeltaJ=double(Q.timeResolved.truthGridDeltaJ);
R.finalTruthInside68=logical(Q.final.truthInside68);
R.timeResolvedTruthInside68=logical(Q.timeResolved.truthInside68);
R.finalTruthInside95=logical(Q.final.truthInside95);
R.timeResolvedTruthInside95=logical(Q.timeResolved.truthInside95);
R.finalL2Profile95Span=finalL2Span;
R.timeResolvedL2Profile95Span=timeL2Span;
R.L2Profile95SpanRatio=spanRatio;
R.I02UnconstrainedAcrossTestedRange=true;
end

% =========================================================================
function assertVectorEqualLocal(actual,expected,tol,label)
a=double(actual(:));
b=double(expected(:));

assert(isequal(size(a),size(b)), ...
    '%s size mismatch.',label);

if isempty(a)
    return
end

assert(max(abs(a-b))<=tol, ...
    '%s numerical mismatch.',label);
end
