function [Mfull,Bbudget] = run_full_model_from_hydraulic_checkpoint_parameterized(checkpointPath,payload)
%RUN_FULL_MODEL_FROM_HYDRAULIC_CHECKPOINT_PARAMETERIZED
% Parameterized sibling of run_full_model_from_hydraulic_checkpoint.
%
% This function preserves the validated reference wrapper architecture while
% externalizing:
%   1) initial soil delta2H by zone
%   2) precipitation delta2H
%   3) leaf-anatomy conditioning
%   4) biochemical policy


if nargin<2 || isempty(payload)
    error('checkpointPath and payload are required.');
end

if ~isfile(checkpointPath)
    error('Hydraulic checkpoint not found: %s',checkpointPath);
end

requireLocal(isstruct(payload),'ParameterizedRunner:BadPayload', ...
    'payload must be a struct.');

requiredPayload = { ...
    'initialSoilDeltaDByZone_permil', ...
    'precipDeltaD_permil', ...
    'leafAnatomyConditioning', ...
    'biochemistryPolicy'};

for k=1:numel(requiredPayload)
    if ~isfield(payload,requiredPayload{k})
        error('Payload is missing required field %s.',requiredPayload{k});
    end
end

Qcheckpoint=load(checkpointPath);
checkpointNames=fieldnames(Qcheckpoint);

for kk=1:numel(checkpointNames)
    eval([checkpointNames{kk} '=Qcheckpoint.(checkpointNames{kk});']); 
end

requiredCheckpointVariables = { ...
    'Crich', ...
    'G76', ...
    'metCoarse', ...
    'plantCoupled', ...
    'timeCoarse', ...
    'xylemTraits', ...
    'leafPV'};

for k = 1:numel(requiredCheckpointVariables)
    if ~exist(requiredCheckpointVariables{k},'var')
        error( ...
            'Checkpoint is missing required variable %s.', ...
            requiredCheckpointVariables{k});
    end
end

timeDays = timeCoarse(:);
nTime = numel(timeDays);
nInterval = nTime - 1;

%% Integrated-run configuration
cfg = struct;
cfg.Rstandard_DH = 155.76e-6;

%% Controlled isotope forcing -- EXTERNALIZED
cfg.forcing.initialSoilDeltaDByZone_permil = ...
    double(payload.initialSoilDeltaDByZone_permil(:));

cfg.forcing.precipDeltaD_permil = ...
    double(payload.precipDeltaD_permil);

cfg.forcing.vaporSource = ...
    'checkpoint metCoarse.deltaDVapor';

%% Soil isotope transport numerical configuration
cfg.soilIsotope = struct;
cfg.soilIsotope.Rstandard_DH = cfg.Rstandard_DH;
cfg.soilIsotope.maxFractionalOutflowPerSubstep = 0.05;
cfg.soilIsotope.minimumSubsteps = 1;
cfg.soilIsotope.maximumSubsteps = 5000;
cfg.soilIsotope.waterTrajectoryTolerance_mol_m2 = 1e-8;
cfg.soilIsotope.rootReleaseMode = 'mixedRootNode';

%% Xylem isotope configuration
cfg.xylemIsotope = xylemTraits;
cfg.xylemIsotope.Rstandard_DH = cfg.Rstandard_DH;
cfg.xylemIsotope.rootIsotopeInterpolation = 'next';
cfg.xylemIsotope.odeRelTol = 1e-8;
cfg.xylemIsotope.odeAbsTol = 1e-10;
cfg.xylemIsotope.odeMaxStepFraction = 1;

%% Spatial leaf isotope traits -- anatomy EXTERNALIZED
cfg.leafIsotope = leafPV;
cfg.leafIsotope.throughflowTolerance_mol_m2_s = 1e-8;
cfg.leafIsotope.epsilonKStomatal = 28.5;
cfg.leafIsotope.epsilonKBoundary = 18.9;

requireLocal(isstruct(payload.leafAnatomyConditioning), ...
    'ParameterizedRunner:BadLeafAnatomy', ...
    'payload.leafAnatomyConditioning must be a struct.');

cfg.leafIsotope.anatomy = payload.leafAnatomyConditioning;

cfg.leafIsotope.outsideXylemPathways.symplasticDiffusiveFraction = 1;
cfg.leafIsotope.isotopeDiffusionConductanceScale = 1;
cfg.leafIsotope.disableIsotopeFractionation = false;

cfg.leafIsotope.spatialSymplast.nCells = 20;
cfg.leafIsotope.spatialSymplast.odeRelTol = 1e-10;
cfg.leafIsotope.spatialSymplast.odeAbsTol = 1e-15;
cfg.leafIsotope.spatialSymplast.odeMaxStepFraction = 0.25;

%% Intracellular water traits
cfg.intracellular = struct;
cfg.intracellular.includePlastid = true;
cfg.intracellular.includeER = false;
cfg.intracellular.waterFractionOfSymplast.cytosol = 0.05;
cfg.intracellular.waterFractionOfSymplast.plastid = 0.003;
cfg.intracellular.exchangeTimescale_s.symplastCytosol = 0.5;
cfg.intracellular.exchangeTimescale_s.cytosolPlastid = 0.1;
cfg.intracellular.sourcePool = 'symplastBulk';
cfg.intracellular.sourceInterpolation = 'linear';
cfg.intracellular.initialState = 'sourceEqual';
cfg.intracellular.odeRelTol = 1e-10;
cfg.intracellular.odeAbsTol = 1e-14;
cfg.intracellular.odeMaxStepFraction = 1;

%% Biochemical diagnostic state -- EXTERNALIZED
requireLocal(isstruct(payload.biochemistryPolicy), ...
    'ParameterizedRunner:BadBiochemistry', ...
    'payload.biochemistryPolicy must be a struct.');

cfg.biochemistry = payload.biochemistryPolicy;

%% 1. Resolved soil isotope transport on saved Richards trajectory
zoneIndex = Crich.soil.zoneIndex(:);
nZone = max(zoneIndex);

if numel(cfg.forcing.initialSoilDeltaDByZone_permil) < nZone
    error( ...
        'Initial soil isotope forcing does not cover all saved soil zones.');
end

deltaDInitialCell_permil = ...
    cfg.forcing.initialSoilDeltaDByZone_permil(zoneIndex);

deltaDPrecipInterval_permil = ...
    expandIntervalForcingLocal( ...
        cfg.forcing.precipDeltaD_permil, ...
        nInterval, ...
        'cfg.forcing.precipDeltaD_permil');

Isoil = ...
    soil_isotope_transport_richards( ...
        Crich.soil, ...
        deltaDInitialCell_permil, ...
        deltaDPrecipInterval_permil, ...
        [], ...
        cfg.soilIsotope);

%% 2. Common-root isotope boundary
Broot = ...
    prepare_root_xylem_isotope_boundary( ...
        Isoil, ...
        Crich, ...
        struct('Rstandard_DH',cfg.Rstandard_DH));

%% 3. Exact-D xylem transport/storage
cfg.xylemIsotope.storageLumenInitialDeltaD = ...
    Broot.initialDeltaD_permil;

cfg.xylemIsotope.storageLivingInitialDeltaD = ...
    Broot.initialDeltaD_permil;

cfg.xylemIsotope.throughflowTolerance_mol_m2_s = 1e-8;

Xxylem = ...
    xylem_transport_dynamic_exactD( ...
        timeDays, ...
        Broot.deltaDNode_permil, ...
        Crich.hydro, ...
        plantCoupled, ...
        cfg.xylemIsotope);

%% 4. Spatial exact-D leaf isotope model
Lleaf = ...
    leaf_hydraulics_dynamic_spatial_symplast( ...
        timeDays, ...
        Xxylem.leafBoundary, ...
        Crich.hydro, ...
        metCoarse, ...
        cfg.leafIsotope);

%% 5. Exact-D intracellular water
Icell = ...
    intracellular_water_dynamic_exactD( ...
        timeDays, ...
        Lleaf, ...
        cfg.intracellular);

%% 6. Biochemical / wax bridge
Bio = ...
    wax_biochemistry_from_intracellular_nominal( ...
        timeDays, ...
        Icell, ...
        cfg.biochemistry);

%% 7. Assemble integrated model object
Mfull = struct;

Mfull.version = ...
    '2026-09-20_full_model_checkpoint_parameterized_v1';

Mfull.timeDays = timeDays;
Mfull.config = cfg;

Mfull.hydraulicCheckpointPath = checkpointPath;
Mfull.hydraulics = Crich;
Mfull.grid = G76;
Mfull.meteorology = metCoarse;

Mfull.soilIsotopes = Isoil;
Mfull.rootBoundary = Broot;
Mfull.xylemIsotopes = Xxylem;
Mfull.leafIsotopes = Lleaf;
Mfull.intracellular = Icell;
Mfull.biochemistry = Bio;

Mfull.scope = { ...
    'precipitation forcing', ...
    'Richards soil isotope transport', ...
    'root-source mixing', ...
    'xylem transport/storage memory', ...
    'leaf evaporative enrichment and spatial Peclet behavior', ...
    'intracellular exchange/storage', ...
    'biochemical source/exchange/incorporation', ...
    'wax homolog isotope composition'};

Mfull.parameterization = struct;
Mfull.parameterization.mode = 'explicitPayload';
Mfull.parameterization.sourceXylemMechanisticInjectionAuthorized = false;
Mfull.parameterization.B2020ExecutionAuthorized = false;

Mfull.knownCurrentOmissions = { ...
    'soil evaporative isotope fractionation is not yet connected', ...
    'accepted-substep Richards internal face isotope fluxes are unavailable in the checkpoint', ...
    'preservation/export isotope effects are not yet implemented', ...
    'B2020 observed xylem delta2H is not yet mapped into the mechanistic soil-root-xylem pathway'};

%% 8. Stage-by-stage isotope/fractionation budget
Bbudget = ...
    build_isotope_fractionation_budget( ...
        Mfull, ...
        deltaDPrecipInterval_permil);

end

function x = expandIntervalForcingLocal(x,nInterval,name)

if isscalar(x)
    x = repmat(x,nInterval,1);
else
    x = x(:);
end

if numel(x) ~= nInterval || any(~isfinite(x))
    error('%s must be scalar or contain exactly nTime-1 finite values.',name);
end

end

function requireLocal(condition,id,message)
if ~condition
    error(id,'%s',message);
end
end
