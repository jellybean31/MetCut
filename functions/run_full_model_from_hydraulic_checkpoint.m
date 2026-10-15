function [Mfull,Bbudget] = run_full_model_from_hydraulic_checkpoint(checkpointPath)
%RUN_FULL_MODEL_FROM_HYDRAULIC_CHECKPOINT
% Parameterized form of driver_full_model using a validated
% hydraulic checkpoint. Isotope, leaf, intracellular, and biochemical
% configuration are otherwise identical to the current production driver.

if nargin<1 || isempty(checkpointPath)
    error('checkpointPath is required.');
end

if ~isfile(checkpointPath)
    error('Hydraulic checkpoint not found: %s',checkpointPath);
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

%% ========================================================================
% Integrated-run configuration
% ========================================================================

cfg = struct;

cfg.Rstandard_DH = 155.76e-6;

% -------------------------------------------------------------------------
% Controlled isotope forcing
% -------------------------------------------------------------------------

cfg.forcing.initialSoilDeltaDByZone_permil = ...
    [-60; -75; -90];

cfg.forcing.precipDeltaD_permil = ...
    -50;

% The saved meteorological vapor forcing remains authoritative for the
% leaf isotope calculation.
cfg.forcing.vaporSource = ...
    'checkpoint metCoarse.deltaDVapor';

% -------------------------------------------------------------------------
% Soil isotope transport numerical configuration
% -------------------------------------------------------------------------

cfg.soilIsotope = struct;
cfg.soilIsotope.Rstandard_DH = cfg.Rstandard_DH;
cfg.soilIsotope.maxFractionalOutflowPerSubstep = 0.05;
cfg.soilIsotope.minimumSubsteps = 1;
cfg.soilIsotope.maximumSubsteps = 5000;
cfg.soilIsotope.waterTrajectoryTolerance_mol_m2 = 1e-8;
cfg.soilIsotope.rootReleaseMode = 'mixedRootNode';

% -------------------------------------------------------------------------
% Xylem isotope configuration
% -------------------------------------------------------------------------

cfg.xylemIsotope = xylemTraits;

cfg.xylemIsotope.Rstandard_DH = cfg.Rstandard_DH;
cfg.xylemIsotope.rootIsotopeInterpolation = 'next';
cfg.xylemIsotope.odeRelTol = 1e-8;
cfg.xylemIsotope.odeAbsTol = 1e-10;
cfg.xylemIsotope.odeMaxStepFraction = 1;


% -------------------------------------------------------------------------
% Spatial leaf isotope traits
% -------------------------------------------------------------------------

cfg.leafIsotope = leafPV;

cfg.leafIsotope.throughflowTolerance_mol_m2_s = 1e-8;

cfg.leafIsotope.epsilonKStomatal = 28.5;
cfg.leafIsotope.epsilonKBoundary = 18.9;

cfg.leafIsotope.anatomy.interveinalDistance_m = 0.0005;
cfg.leafIsotope.anatomy.mesophyllThickness_m = 0.00015;
cfg.leafIsotope.anatomy.effectiveFlowAreaFraction = 0.02;
cfg.leafIsotope.anatomy.flowFraction = 1;
cfg.leafIsotope.anatomy.tortuosity = 1.5;
cfg.leafIsotope.anatomy.lateralPathFraction = 0.5;
cfg.leafIsotope.anatomy.verticalPathFraction = 0.5;

cfg.leafIsotope.outsideXylemPathways.symplasticDiffusiveFraction = 1;
cfg.leafIsotope.isotopeDiffusionConductanceScale = 1;
cfg.leafIsotope.disableIsotopeFractionation = false;

cfg.leafIsotope.spatialSymplast.nCells = 20;
cfg.leafIsotope.spatialSymplast.odeRelTol = 1e-10;
cfg.leafIsotope.spatialSymplast.odeAbsTol = 1e-15;
cfg.leafIsotope.spatialSymplast.odeMaxStepFraction = 0.25;

% -------------------------------------------------------------------------
% Intracellular water traits
% -------------------------------------------------------------------------

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

% -------------------------------------------------------------------------
% Biochemical diagnostic state
% -------------------------------------------------------------------------

cfg.biochemistry = struct;
cfg.biochemistry.parameterMode = 'coreMidpointDiagnostic';
cfg.biochemistry.routeFractions = struct('C16',0,'C18',1);
cfg.biochemistry.alkaneCarbonNumber = [27 29 31 33];

%% ========================================================================
% 1. Resolved soil isotope transport on saved Richards trajectory
% ========================================================================

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

%% ========================================================================
% 2. Common-root isotope boundary
% ========================================================================

Broot = ...
    prepare_root_xylem_isotope_boundary( ...
        Isoil, ...
        Crich, ...
        struct('Rstandard_DH',cfg.Rstandard_DH));

%% ========================================================================
% 3. Exact-D xylem transport/storage
% ========================================================================

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

%% ========================================================================
% 4. Spatial exact-D leaf isotope model
% ========================================================================

Lleaf = ...
    leaf_hydraulics_dynamic_spatial_symplast( ...
        timeDays, ...
        Xxylem.leafBoundary, ...
        Crich.hydro, ...
        metCoarse, ...
        cfg.leafIsotope);

%% ========================================================================
% 5. Exact-D intracellular water
% ========================================================================

Icell = ...
    intracellular_water_dynamic_exactD( ...
        timeDays, ...
        Lleaf, ...
        cfg.intracellular);

%% ========================================================================
% 6. Biochemical / wax bridge
% ========================================================================

Bio = ...
    wax_biochemistry_from_intracellular_nominal( ...
        timeDays, ...
        Icell, ...
        cfg.biochemistry);

%% ========================================================================
% 7. Assemble integrated model object
% ========================================================================

Mfull = struct;

Mfull.version = ...
    '2026-08-18_full_model_checkpoint_integration_v1';

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

Mfull.knownCurrentOmissions = { ...
    'soil evaporative isotope fractionation is not yet connected', ...
    'substep Richards internal face isotope fluxes are unavailable', ...
    'biochemical trajectory uses midpoint', ...
    'preservation/export isotope effects are not yet implemented'};

%% ========================================================================
% 8. Stage-by-stage isotope/fractionation budget
% ========================================================================

Bbudget = ...
    build_isotope_fractionation_budget( ...
        Mfull, ...
        deltaDPrecipInterval_permil);

%% ========================================================================
% Console checkpoint
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('FULL MODEL CHECKPOINT INTEGRATION V1\n');
fprintf('============================================================\n');
fprintf('Hydraulics: 76-cell checkpoint\n');
fprintf('Timesteps:  %d\n',nTime);
fprintf('Soil cells: %d\n',size(Isoil.deltaD_permil,1));
fprintf('Leaf spatial cells: %d\n',cfg.leafIsotope.spatialSymplast.nCells);
fprintf('Biochemical state: %s\n',cfg.biochemistry.parameterMode);
fprintf('------------------------------------------------------------\n');
fprintf('Maximum soil D closure:          %.12g mol D m^-2\n', ...
    max(abs(Isoil.stepIsotopeBalanceResidual_mol_D_m2),[],'omitnan'));
fprintf('Maximum xylem D closure:         %.12g mol D m^-2\n', ...
    max(abs(Xxylem.deuteriumInventoryClosureResidual_molD_m2),[],'omitnan'));
fprintf('Maximum leaf D closure:          %.12g mol D m^-2\n', ...
    max(abs(Lleaf.DInventoryClosureResidual_molD_m2),[],'omitnan'));
fprintf('Maximum intracellular D closure: %.12g mol D m^-2\n', ...
    max(abs(Icell.wholeSystemDClosureResidual_molD_m2),[],'omitnan'));
fprintf('------------------------------------------------------------\n');
fprintf('Final root-source deltaD:        %.6f per mil\n', ...
    Broot.deltaDNode_permil(end));
fprintf('Final xylem-delivered deltaD:    %.6f per mil\n', ...
    Xxylem.leafBoundary.deltaDWater(end));
fprintf('Final bulk symplast deltaD:      %.6f per mil\n', ...
    Lleaf.deltaDSymplastBulk(end));
fprintf('Final cytosol deltaD:            %.6f per mil\n', ...
    Icell.deltaD.cytosol(end));
fprintf('Final plastid deltaD:            %.6f per mil\n', ...
    Icell.deltaD.plastid(end));
fprintf('Final C29 wax deltaD:            %.6f per mil\n', ...
    Bio.wax.deltaDMatrix(2,end));
fprintf('============================================================\n');

disp(Bbudget.transitionMetadata);

end

%% ========================================================================
% Local helper
% ========================================================================

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
