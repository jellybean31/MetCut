function R = test_METCUT_hydraulic_solver_regression()
%TEST_METCUT_HYDRAULIC_SOLVER_REGRESSION
%
% METCUT integration regression for the coupled hydraulic
% solver. The test reruns the first 24 h (25 hourly nodes) of the
% 76-cell reference trajectory using the package-local hydraulic checkpoint
% as the fixed-point initial guess and numerical reference.
%
% This exercises:
%   couple_richards_soil_plant_hydraulics
%   solve_plant_hydraulics
%   soil_water_balance_richards_1d_prescribed_root


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
assert(isfile(hydPath),'Package-local hydraulic reference missing.');

Q=load(hydPath);

requiredCheckpointFields={ ...
    'Crich','G76','metCoarse','plantCoupled','precipCoarse', ...
    'richardsCfg','root76','soil76','timeCoarse','xylemTraits','leafPV', ...
    'couplingCfgRichards'};

for i=1:numel(requiredCheckpointFields)
    assert(isfield(Q,requiredCheckpointFields{i}), ...
        'Hydraulic reference lacks required field %s.',requiredCheckpointFields{i});
end

assert(Q.G76.nCell==76,'Expected accepted 76-cell hydraulic grid.');

% -------------------------------------------------------------------------
% Short execution window: first 24 h / 25 nodes.
% -------------------------------------------------------------------------
nNode=25;
assert(numel(Q.timeCoarse)>=nNode, ...
    'Hydraulic reference is shorter than the requested regression window.');

idx=1:nNode;
timeDays=double(Q.timeCoarse(idx));
psiRef=double(Q.Crich.psiSoilCells_MPa(:,idx));

assert(isequal(size(psiRef),[Q.G76.nCell nNode]), ...
    'Accepted soil-pressure reference has unexpected dimensions.');

psiInitialCell=psiRef(:,1);
psiInitialGuess=psiRef;

precip=subsetNodeVectorLocal(Q.precipCoarse,idx,numel(Q.timeCoarse));
met=subsetMetLocal(Q.metCoarse,idx,numel(Q.timeCoarse));

tRun=tic;
C=couple_richards_soil_plant_hydraulics( ...
    timeDays, ...
    psiInitialCell, ...
    psiInitialGuess, ...
    precip, ...
    met, ...
    Q.root76, ...
    Q.plantCoupled, ...
    Q.xylemTraits, ...
    Q.leafPV, ...
    Q.soil76, ...
    Q.G76, ...
    Q.richardsCfg, ...
    Q.couplingCfgRichards, ...
    struct);
runtime_s=toc(tRun);

% -------------------------------------------------------------------------
% Coupled fixed-point / execution gates.
% -------------------------------------------------------------------------
assert(C.converged,'Short coupled hydraulic regression did not converge.');
assert(C.nIterations>=2, ...
    'Expected at least two coupling iterations because first dJ is Inf.');
assert(C.nAttempts>=C.nIterations,'Hydraulic attempt bookkeeping regression.');
assert(C.nSoilCells==76,'Hydraulic soil-cell-count regression.');

if isfield(Q.couplingCfgRichards,'psiTolerance_MPa') && ...
        ~isempty(Q.couplingCfgRichards.psiTolerance_MPa)
    psiTolerance_MPa=Q.couplingCfgRichards.psiTolerance_MPa;
else
    psiTolerance_MPa=1e-6;
end

if isfield(Q.couplingCfgRichards,'areaBasisTolerance_mol_m2_s') && ...
        ~isempty(Q.couplingCfgRichards.areaBasisTolerance_mol_m2_s)
    areaBasisTolerance_mol_m2_s= ...
        Q.couplingCfgRichards.areaBasisTolerance_mol_m2_s;
else
    areaBasisTolerance_mol_m2_s=1e-12;
end

assert(C.maxFinalPsiResidual_MPa<=psiTolerance_MPa, ...
    'Final coupled soil-pressure residual exceeds configured/default tolerance.');
assert(C.maxAreaBasisResidual_mol_m2_s<=areaBasisTolerance_mol_m2_s, ...
    'Leaf-area/ground-area root-flux conversion regression.');

H=C.hydro;
S=C.soil;

assert(~any(H.solverFallback), ...
    'Plant hydraulic solver used fallback within the 24-h window.');
assert(all(H.hydraulicallyFeasible), ...
    'Plant hydraulic solver reported an infeasible state.');

% -------------------------------------------------------------------------
% Richards mass conservation.
% -------------------------------------------------------------------------
soilStep=max(abs(double(S.stepMassBalanceResidual_mol_m2)),[],'omitnan');
soilCum=max(abs(double(S.cumulativeMassBalanceResidual_mol_m2)),[],'omitnan');

assert(soilStep<=1e-9,'Richards step mass-balance regression.');
assert(soilCum<=1e-9,'Richards cumulative mass-balance regression.');

% -------------------------------------------------------------------------
% Plant hydraulic closure.
% -------------------------------------------------------------------------
wholePlantClosure=max(abs(double( ...
    H.wholePlantStorageFluxClosureResidual)),[],'omitnan');
xylemStorageClosure=max(abs(double( ...
    H.storageFluxClosureResidual)),[],'omitnan');
leafStorageClosure=max(abs(double( ...
    H.leafStorageFluxClosureResidual)),[],'omitnan');
stemFlowClosure=max(abs(double(H.stemFlowResidual)),[],'omitnan');

assert(wholePlantClosure<=1e-8, ...
    'Whole-plant hydraulic storage/flow closure regression.');
assert(xylemStorageClosure<=1e-8, ...
    'Xylem storage-flow closure regression.');
assert(leafStorageClosure<=1e-8, ...
    'Leaf storage-flow closure regression.');
assert(stemFlowClosure<=1e-8, ...
    'Stem-flow closure regression.');

% -------------------------------------------------------------------------
% Accepted trajectory regression.
%
% Because the reference soil trajectory is supplied only as the initial
% fixed-point guess, agreement here still requires a complete plant ->
% root-flux -> Richards -> fixed-point evaluation by the current code.
% -------------------------------------------------------------------------
soilPsiResidual=max(abs(double(C.psiSoilCells_MPa)-psiRef),[],'all','omitnan');

Href=Q.Crich.hydro;

psiRootRef=double(Href.psiRoot(idx));
psiXylemRef=double(Href.psiXylem(idx));
psiLeafRef=double(Href.psiLeaf(idx));

psiRootResidual=max(abs(double(H.psiRoot)-psiRootRef),[],'omitnan');
psiXylemResidual=max(abs(double(H.psiXylem)-psiXylemRef),[],'omitnan');
psiLeafResidual=max(abs(double(H.psiLeaf)-psiLeafRef),[],'omitnan');

% Coupled reference was itself accepted at 1e-6 MPa fixed-point tolerance.
% Use the same scale for release trajectory comparison.
trajectoryTolerance_MPa=2e-6;

assert(soilPsiResidual<=trajectoryTolerance_MPa, ...
    '24-h soil-pressure trajectory regression.');
assert(psiRootResidual<=trajectoryTolerance_MPa, ...
    '24-h root-pressure trajectory regression.');
assert(psiXylemResidual<=trajectoryTolerance_MPa, ...
    '24-h xylem-pressure trajectory regression.');
assert(psiLeafResidual<=trajectoryTolerance_MPa, ...
    '24-h leaf-pressure trajectory regression.');

% -------------------------------------------------------------------------
% Physical state gates.
% -------------------------------------------------------------------------
assert(all(isfinite(C.psiSoilCells_MPa),'all'), ...
    'Nonfinite coupled soil-pressure state.');
assert(all(isfinite(H.psiRoot)), ...
    'Nonfinite root pressure.');
assert(all(isfinite(H.psiXylem)), ...
    'Nonfinite xylem pressure.');
assert(all(isfinite(H.psiLeaf)), ...
    'Nonfinite leaf pressure.');
assert(all(isfinite(H.E)) && all(H.E>=0), ...
    'Invalid transpiration trajectory.');
assert(all(isfinite(H.gs)) && all(H.gs>=0), ...
    'Invalid stomatal-conductance trajectory.');

if H.xylemCapacitanceEnabled
    A=double(H.hydraulicStorageAmount_mol_m2);
    assert(all(isfinite(A)) && all(A>0), ...
        'Invalid xylem hydraulic-storage trajectory.');
    assert(~any(H.hydraulicStorageCapacityViolated), ...
        'Xylem hydraulic-storage capacity violation.');
end

if H.leafCapacitanceEnabled
    A=double(H.leafHydraulicStorageAmount_mol_m2);
    assert(all(isfinite(A)) && all(A>0), ...
        'Invalid leaf hydraulic-storage trajectory.');
end

% -------------------------------------------------------------------------
% Result record.
% -------------------------------------------------------------------------
R=struct;
R.pass=true;
R.runtime_s=runtime_s;
R.nNodes=nNode;
R.nSoilCells=C.nSoilCells;
R.nIterations=C.nIterations;
R.nAttempts=C.nAttempts;
R.finalRelaxation=C.finalRelaxation;
R.maxFinalPsiResidual_MPa=C.maxFinalPsiResidual_MPa;
R.maxAreaBasisResidual_mol_m2_s=C.maxAreaBasisResidual_mol_m2_s;
R.psiTolerance_MPa=psiTolerance_MPa;
R.areaBasisTolerance_mol_m2_s=areaBasisTolerance_mol_m2_s;

R.richardsStepMassResidual_mol_m2=soilStep;
R.richardsCumulativeMassResidual_mol_m2=soilCum;
R.wholePlantClosure_mol_m2_s=wholePlantClosure;
R.xylemStorageClosure_mol_m2_s=xylemStorageClosure;
R.leafStorageClosure_mol_m2_s=leafStorageClosure;
R.stemFlowClosure_mol_m2_s=stemFlowClosure;

R.soilPsiReferenceResidual_MPa=soilPsiResidual;
R.psiRootReferenceResidual_MPa=psiRootResidual;
R.psiXylemReferenceResidual_MPa=psiXylemResidual;
R.psiLeafReferenceResidual_MPa=psiLeafResidual;
R.trajectoryTolerance_MPa=trajectoryTolerance_MPa;
end

% =========================================================================
function x=subsetNodeVectorLocal(x,idx,nFull)
x0=x;
if isscalar(x0)
    x=x0;
    return
end

x0=x0(:).';

if numel(x0)==nFull
    x=x0(idx);
elseif numel(x0)==nFull-1
    % Interval forcing. For N selected nodes retain N-1 intervals.
    x=x0(idx(1:end-1));
else
    error('Unexpected forcing-vector length in hydraulic reference.');
end
end

% =========================================================================
function M=subsetMetLocal(M0,idx,nFull)
M=struct;
names=fieldnames(M0);

for i=1:numel(names)
    nm=names{i};
    x=M0.(nm);

    if ~isnumeric(x) && ~islogical(x)
        M.(nm)=x;
        continue
    end

    if isscalar(x)
        M.(nm)=x;
        continue
    end

    sz=size(x);

    if isvector(x) && numel(x)==nFull
        wasColumn=iscolumn(x);
        y=x(:).';
        y=y(idx);
        if wasColumn, y=y(:); end
        M.(nm)=y;
    elseif sz(end)==nFull
        subs=repmat({':'},1,ndims(x));
        subs{end}=idx;
        M.(nm)=x(subs{:});
    elseif sz(1)==nFull
        subs=repmat({':'},1,ndims(x));
        subs{1}=idx;
        M.(nm)=x(subs{:});
    else
        M.(nm)=x;
    end
end
end
