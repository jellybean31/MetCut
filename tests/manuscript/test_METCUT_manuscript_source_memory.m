function R = test_METCUT_manuscript_source_memory()
%TEST_METCUT_MANUSCRIPT_SOURCE_MEMORY
%
% Regression for:
%   1. explicit isotope-source override replay;
%   2. fixed-hydraulic temporal-memory replay.
%
% The first gate is an equivalence test:
% source tagging at baseline source values must reproduce the ordinary
% isotope-forcing + wax-timing replay.
%
% The second gate perturbs a single source class and verifies that isotope
% outputs respond.
%
% The third gate evaluates two temporal-memory settings and verifies
% that the I01/I02 memory controls are wired into isotope propagation.

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

[M0,~]=run_full_model_from_hydraulic_checkpoint(hydPath);

nInt=numel(M0.timeDays)-1;
d0=expandIntervalLocal( ...
    M0.config.forcing.precipDeltaD_permil,nInt);

source0=struct;
source0.precipDeltaDInterval_permil=d0;
source0.initialSoilDeltaDByZone_permil= ...
    double(M0.config.forcing.initialSoilDeltaDByZone_permil(:));

assert(isfield(M0,'rootBoundary') && ...
    isfield(M0.rootBoundary,'initialDeltaD_permil') && ...
    isfinite(M0.rootBoundary.initialDeltaD_permil), ...
    'Baseline root-source initial deltaD is unavailable.');

dStorage0=double(M0.rootBoundary.initialDeltaD_permil);
source0.xylemStorageLumenInitialDeltaD_permil=dStorage0;
source0.xylemStorageLivingInitialDeltaD_permil=dStorage0;

% -------------------------------------------------------------------------
% 1. Baseline source-override equivalence.
% -------------------------------------------------------------------------
[Mstd,Bstd]=rerun_full_model_isotope_forcing_with_wax_timing( ...
    M0,d0,[],[]);

[Msrc,Bsrc]=rerun_full_model_isotope_source_overrides_with_wax_timing( ...
    M0,source0,[],[]);

compareLocal(Mstd.soilIsotopes.deltaD_permil, ...
    Msrc.soilIsotopes.deltaD_permil,'soil isotope state');
compareLocal(Mstd.rootBoundary.deltaDNode_permil, ...
    Msrc.rootBoundary.deltaDNode_permil,'root-source state');
compareLocal(Mstd.xylemIsotopes.leafBoundary.deltaDWater, ...
    Msrc.xylemIsotopes.leafBoundary.deltaDWater,'xylem delivered state');
compareLocal(Mstd.leafIsotopes.deltaDSymplastBulk, ...
    Msrc.leafIsotopes.deltaDSymplastBulk,'bulk symplast state');
compareLocal(Mstd.intracellular.deltaD.cytosol, ...
    Msrc.intracellular.deltaD.cytosol,'cytosol state');
compareLocal(Mstd.intracellular.deltaD.plastid, ...
    Msrc.intracellular.deltaD.plastid,'plastid state');
compareLocal(Mstd.biochemistry.waxStageExpanded.deltaDMatrix, ...
    Msrc.biochemistry.waxStageExpanded.deltaDMatrix, ...
    'architecture-v2 instantaneous plant wax');

assert(isequaln(Msrc.hydraulics,M0.hydraulics), ...
    'Source-override replay modified hydraulics.');

assert(isfield(Msrc,'sourceTagging') && ...
    string(Msrc.sourceTagging.version)=="2026-08-30_H8_source_tagging_v1", ...
    'Source-tagging metadata regression.');

% -------------------------------------------------------------------------
% 2. Single-source perturbation: pre-existing xylem lumen storage.
% -------------------------------------------------------------------------
sourceLumen=source0;
sourceLumen.xylemStorageLumenInitialDeltaD_permil= ...
    source0.xylemStorageLumenInitialDeltaD_permil+100;

[Mlumen,Blumen]=rerun_full_model_isotope_source_overrides_with_wax_timing( ...
    M0,sourceLumen,[],[]);

assert(isequaln(Mlumen.hydraulics,M0.hydraulics), ...
    'Source perturbation modified hydraulics.');

dX=max(abs(double( ...
    Mlumen.xylemIsotopes.leafBoundary.deltaDWater(:))-double( ...
    Msrc.xylemIsotopes.leafBoundary.deltaDWater(:))),[],'omitnan');

dWax=max(abs(double( ...
    Mlumen.biochemistry.waxStageExpanded.deltaDMatrix(:))-double( ...
    Msrc.biochemistry.waxStageExpanded.deltaDMatrix(:))),[],'omitnan');

assert(isfinite(dX) && dX>1e-8, ...
    'Lumen-source perturbation did not propagate to xylem output.');
assert(isfinite(dWax) && dWax>1e-8, ...
    'Lumen-source perturbation did not propagate to plant-wax output.');

checkClosureLocal(Bstd,'standard baseline');
checkClosureLocal(Bsrc,'source-tag baseline');
checkClosureLocal(Blumen,'lumen-source perturbation');

% -------------------------------------------------------------------------
% 3. Temporal-memory wiring.
%
% Accepted synthetic truth used by the full-chain wax-only recovery:
%   I01 mobileTransitDays       = 2.03690139 d
%   I02 storageExchangeTauDays = 5.38306961 d
%
% The comparison case uses the accepted replayable upper bounds:
%   I01 = 5 d
%   I02 = 20 d
% -------------------------------------------------------------------------
memoryTruth=struct( ...
    'mobileTransitDays',2.03690139, ...
    'storageExchangeTauDays',5.38306961);

memoryLong=struct( ...
    'mobileTransitDays',5, ...
    'storageExchangeTauDays',20);

[MmemTruth,BmemTruth]=rerun_full_model_temporal_memory_source_case( ...
    M0,source0,memoryTruth,[],[]);

[MmemLong,BmemLong]=rerun_full_model_temporal_memory_source_case( ...
    M0,source0,memoryLong,[],[]);

assert(logical(MmemTruth.temporalMemoryCase.hydraulicsReused), ...
    'Temporal-memory truth case did not report hydraulic reuse.');
assert(logical(MmemLong.temporalMemoryCase.hydraulicsReused), ...
    'Temporal-memory long case did not report hydraulic reuse.');

assert(isequaln(MmemTruth.hydraulics,M0.hydraulics) && ...
       isequaln(MmemLong.hydraulics,M0.hydraulics), ...
    'Temporal-memory replay modified hydraulics.');

assert(abs(MmemTruth.temporalMemoryCase.mobileTransitDays- ...
    memoryTruth.mobileTransitDays)<=1e-12, ...
    'I01 temporal-memory metadata regression.');
assert(abs(MmemTruth.temporalMemoryCase.storageExchangeTauDays- ...
    memoryTruth.storageExchangeTauDays)<=1e-12, ...
    'I02 temporal-memory metadata regression.');

dMemoryX=max(abs(double( ...
    MmemLong.xylemIsotopes.leafBoundary.deltaDWater(:))-double( ...
    MmemTruth.xylemIsotopes.leafBoundary.deltaDWater(:))),[],'omitnan');

dMemoryWax=max(abs(double( ...
    MmemLong.biochemistry.waxStageExpanded.deltaDMatrix(:))-double( ...
    MmemTruth.biochemistry.waxStageExpanded.deltaDMatrix(:))),[],'omitnan');

assert(isfinite(dMemoryX) && dMemoryX>1e-8, ...
    'Temporal-memory perturbation did not alter xylem isotope output.');
assert(isfinite(dMemoryWax) && dMemoryWax>1e-8, ...
    'Temporal-memory perturbation did not alter plant-wax output.');

checkClosureLocal(BmemTruth,'temporal-memory truth');
checkClosureLocal(BmemLong,'temporal-memory long');

% Record C29 contrasts for later numerical fingerprinting.
[Ctruth,Dtruth]=waxMatrixLocal(MmemTruth);
[Clong,Dlong]=waxMatrixLocal(MmemLong);
assert(isequal(Ctruth,Clong),'Temporal-memory homolog-set regression.');

i29=find(Ctruth==29,1);
assert(~isempty(i29),'C29 missing from temporal-memory output.');

finalC29Truth=Dtruth(i29,end);
finalC29Long=Dlong(i29,end);
finalC29Difference=finalC29Long-finalC29Truth;

fprintf('\n');
fprintf('============================================================\n');
fprintf('METCUT SOURCE / TEMPORAL-MEMORY REGRESSION\n');
fprintf('============================================================\n');
fprintf('Baseline source-override equivalence: PASS\n');
fprintf('Lumen +100 permil -> max xylem response: %.12g permil\n',dX);
fprintf('Lumen +100 permil -> max wax response:   %.12g permil\n',dWax);
fprintf('Memory truth I01/I02:                    %.8f / %.8f d\n', ...
    memoryTruth.mobileTransitDays,memoryTruth.storageExchangeTauDays);
fprintf('Memory long I01/I02:                     %.8f / %.8f d\n', ...
    memoryLong.mobileTransitDays,memoryLong.storageExchangeTauDays);
fprintf('Memory contrast max xylem response:      %.12g permil\n',dMemoryX);
fprintf('Memory contrast max wax response:        %.12g permil\n',dMemoryWax);
fprintf('Final C29 truth / long / difference:     %.12g / %.12g / %.12g permil\n', ...
    finalC29Truth,finalC29Long,finalC29Difference);
fprintf('PASS\n');
fprintf('============================================================\n');

R=struct;
R.pass=true;
R.baselineSourceOverrideEquivalent=true;
R.baselineStorageDeltaD_permil=dStorage0;
R.lumenPerturbation_permil=100;
R.lumenPerturbationMaxXylemResponse_permil=dX;
R.lumenPerturbationMaxWaxResponse_permil=dWax;
R.memoryTruth=memoryTruth;
R.memoryLong=memoryLong;
R.memoryMaxXylemDifference_permil=dMemoryX;
R.memoryMaxWaxDifference_permil=dMemoryWax;
R.finalC29Truth_permil=finalC29Truth;
R.finalC29Long_permil=finalC29Long;
R.finalC29Difference_permil=finalC29Difference;
end

% =========================================================================
function x=expandIntervalLocal(x,nInt)
if isscalar(x)
    x=repmat(double(x),nInt,1);
else
    x=double(x(:));
end
assert(numel(x)==nInt && all(isfinite(x)), ...
    'Baseline precipitation forcing must resolve to nTime-1 finite values.');
end

% =========================================================================
function compareLocal(a,b,label)
a=double(a);
b=double(b);
assert(isequal(size(a),size(b)),'%s size regression.',label);
assert(isequal(isnan(a),isnan(b)),'%s NaN-pattern regression.',label);
assert(isequal(isinf(a),isinf(b)),'%s Inf-pattern regression.',label);

m=isfinite(a)&isfinite(b);
if any(m(:))
    d=abs(a(m)-b(m));
    scale=max(abs(a(m)),abs(b(m)));
    rel=d./max(scale,1);
    assert(max(d)<=1e-10 && max(rel)<=1e-12, ...
        '%s source-override equivalence regression.',label);
end
end

% =========================================================================
function checkClosureLocal(B,label)
assert(B.closure.maxSoilDResidual_molD_m2<1e-10, ...
    '%s soil exact-D closure regression.',label);
assert(B.closure.maxXylemDResidual_molD_m2<1e-10, ...
    '%s xylem exact-D closure regression.',label);
assert(B.closure.maxLeafDResidual_molD_m2<1e-9, ...
    '%s leaf exact-D closure regression.',label);
assert(B.closure.maxIntracellularDResidual_molD_m2<1e-12, ...
    '%s intracellular exact-D closure regression.',label);
end

% =========================================================================
function [C,D]=waxMatrixLocal(M)
W=M.biochemistry.waxStageExpanded;
C=double(W.alkaneCarbonNumber(:));
D=double(W.deltaDMatrix);
if size(D,1)~=numel(C) && size(D,2)==numel(C)
    D=D.';
end
assert(size(D,1)==numel(C),'Wax matrix/carbon-number shape regression.');
end
