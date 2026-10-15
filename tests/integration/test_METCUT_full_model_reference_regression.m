function R = test_METCUT_full_model_reference_regression()
%TEST_METCUT_FULL_MODEL_REFERENCE_REGRESSION
%
% METCUT integration regression for the instantaneous forward model.


thisFile=mfilename('fullpath');
integrationDir=fileparts(thisFile);
testsDir=fileparts(integrationDir);
rootDir=fileparts(testsDir);
fnDir=fullfile(rootDir,'functions');
refDir=fullfile(rootDir,'data','reference');

originalPath=path;
cleanup=onCleanup(@() path(originalPath)); 
addpath(fnDir,'-begin');
rehash

hydPath=fullfile(refDir, ...
    'hydraulic_checkpoint_76cell_8day_hourly_corrected_mass_conservative_storagefixed.mat');
valPath=fullfile(refDir,'model_validation_acceptance_full_v3.mat');

assert(isfile(hydPath),'Package-local hydraulic reference missing.');
assert(isfile(valPath),'Package-local validation fingerprint missing.');

[M,B]=run_full_model_from_hydraulic_checkpoint(hydPath);

soilStep=max(abs(double( ...
    M.soilIsotopes.stepIsotopeBalanceResidual_mol_D_m2(:))),[],'omitnan');
soilCum=max(abs(double( ...
    M.soilIsotopes.cumulativeIsotopeBalanceResidual_mol_D_m2(:))),[],'omitnan');
xylemD=max(abs(double( ...
    M.xylemIsotopes.deuteriumInventoryClosureResidual_molD_m2(:))),[],'omitnan');
leafD=max(abs(double( ...
    M.leafIsotopes.DInventoryClosureResidual_molD_m2(:))),[],'omitnan');
intracellularD=max(abs(double( ...
    M.intracellular.wholeSystemDClosureResidual_molD_m2(:))),[],'omitnan');

assert(soilStep<1e-10,'Soil step exact-D closure regression.');
assert(soilCum<1e-10,'Soil cumulative exact-D closure regression.');
assert(xylemD<1e-10,'Xylem exact-D closure regression.');
assert(leafD<1e-9,'Leaf exact-D closure regression.');
assert(intracellularD<1e-12,'Intracellular exact-D closure regression.');

W=double(M.biochemistry.wax.deltaDMatrix);

% The accepted baseline Bio.wax schema does not necessarily carry explicit
% carbon numbers. When absent, the production baseline row order is the
% validated C27/C29/C31/C33 order used by the historical acceptance suite.
if isfield(M.biochemistry.wax,'alkaneCarbonNumber')
    C=double(M.biochemistry.wax.alkaneCarbonNumber(:).');
else
    assert(size(W,1)==4 || size(W,2)==4, ...
        'Wax carbon numbers absent and four-homologue fallback is invalid.');
    C=[27 29 31 33];
end

if size(W,1)~=numel(C) && size(W,2)==numel(C)
    W=W.';
end

assert(size(W,1)==numel(C), ...
    'Wax matrix/carbon-number shape regression.');
assert(all(isfinite(W(:))),'Nonfinite instantaneous plant-wax output.');

finalRoot=lastFiniteLocal(M.rootBoundary.deltaDNode_permil);
finalXylem=lastFiniteLocal(M.xylemIsotopes.leafBoundary.deltaDWater);
finalLeaf=lastFiniteLocal(M.leafIsotopes.deltaDSymplastBulk);
finalCy=lastFiniteLocal(M.intracellular.deltaD.cytosol);
finalPl=lastFiniteLocal(M.intracellular.deltaD.plastid);

finalWax=NaN(1,4);
targets=[27 29 31 33];
for i=1:4
    idx=find(C==targets(i),1);
    assert(~isempty(idx),'Missing expected alkane homologue C%d.',targets(i));
    finalWax(i)=lastFiniteLocal(W(idx,:));
end

% Accepted package-local release fingerprints from A05b-v4.
tolState=1e-6;
tolWax=1e-6;

assert(abs(finalRoot - (-73.2800577903723))<=tolState,'Final root dD regression.');
assert(abs(finalXylem - (-66.8297768879364))<=tolState,'Final xylem dD regression.');
assert(abs(finalLeaf - (-25.6304872328379))<=tolState,'Final leaf dD regression.');
assert(abs(finalCy - (-25.6304871371434))<=tolState,'Final cytosol dD regression.');
assert(abs(finalPl - (-25.6304871190879))<=tolState,'Final plastid dD regression.');

expectedWax=[ ...
    -146.557419798513 ...
    -146.344758391238 ...
    -146.158679648317 ...
    -145.994492513254];

assert(max(abs(finalWax-expectedWax))<=tolWax, ...
    'Final C27/C29/C31/C33 instantaneous wax regression.');

% Historical accepted C29 fingerprint remains within project tiny-regression
% tolerance.
S=load(valPath,'R');
oldC29=double(S.R.outputs{2}.finalC29_permil);
assert(abs(finalWax(2)-oldC29)<=1e-6, ...
    'Historical final-C29 regression.');

% Parameterized sibling must remain identical under exact reference defaults.
P=struct;
P.initialSoilDeltaDByZone_permil=[-60;-75;-90];
P.precipDeltaD_permil=-50;
P.leafAnatomyConditioning=struct( ...
    'interveinalDistance_m',0.0005, ...
    'mesophyllThickness_m',0.00015, ...
    'effectiveFlowAreaFraction',0.02, ...
    'flowFraction',1, ...
    'tortuosity',1.5, ...
    'lateralPathFraction',0.5, ...
    'verticalPathFraction',0.5);
P.biochemistryPolicy=struct;
P.biochemistryPolicy.parameterMode='coreMidpointDiagnostic';
P.biochemistryPolicy.routeFractions=struct('C16',0,'C18',1);
P.biochemistryPolicy.alkaneCarbonNumber=[27 29 31 33];

[Mp,Bp]=run_full_model_from_hydraulic_checkpoint_parameterized(hydPath,P);

compareLocal(M.soilIsotopes.deltaD_permil,Mp.soilIsotopes.deltaD_permil,'soil');
compareLocal(M.rootBoundary.deltaDNode_permil,Mp.rootBoundary.deltaDNode_permil,'root');
compareLocal(M.xylemIsotopes.leafBoundary.deltaDWater, ...
    Mp.xylemIsotopes.leafBoundary.deltaDWater,'xylem');
compareLocal(M.leafIsotopes.deltaDSymplastBulk, ...
    Mp.leafIsotopes.deltaDSymplastBulk,'leaf');
compareLocal(M.intracellular.deltaD.cytosol, ...
    Mp.intracellular.deltaD.cytosol,'cytosol');
compareLocal(M.intracellular.deltaD.plastid, ...
    Mp.intracellular.deltaD.plastid,'plastid');
compareLocal(M.biochemistry.wax.deltaDMatrix, ...
    Mp.biochemistry.wax.deltaDMatrix,'wax');

assert(isequaln(B.transitionMetadata,Bp.transitionMetadata), ...
    'Parameterized/reference budget metadata regression.');

R=struct;
R.pass=true;
R.finalRootDeltaD_permil=finalRoot;
R.finalXylemDeltaD_permil=finalXylem;
R.finalLeafBulkDeltaD_permil=finalLeaf;
R.finalCytosolDeltaD_permil=finalCy;
R.finalPlastidDeltaD_permil=finalPl;
R.finalWaxDeltaD_permil=finalWax;
R.closure=struct('soilStep',soilStep,'soilCumulative',soilCum, ...
    'xylem',xylemD,'leaf',leafD,'intracellular',intracellularD);
end

function compareLocal(a,b,label)
assert(isequal(size(a),size(b)),'%s size regression.',label);
a=double(a); b=double(b);
assert(isequal(isnan(a),isnan(b)),'%s NaN-pattern regression.',label);
assert(isequal(isinf(a),isinf(b)),'%s Inf-pattern regression.',label);
m=isfinite(a)&isfinite(b);
if any(m(:))
    d=abs(a(m)-b(m));
    scale=max(abs(a(m)),abs(b(m)));
    rel=d./max(scale,1);
    assert(max(d)<=1e-10 && max(rel)<=1e-12, ...
        '%s parameterized/reference regression.',label);
end
end

function x=lastFiniteLocal(v)
v=double(v(:));
idx=find(isfinite(v),1,'last');
assert(~isempty(idx),'Expected finite state.');
x=v(idx);
end
