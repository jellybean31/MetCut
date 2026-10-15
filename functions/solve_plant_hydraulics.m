function hydro = solve_plant_hydraulics( ...
    timeDays, soil, met, rootTraits, plantState, ...
    xylemTraits, leafTraits, initialState)
%SOLVE_PLANT_HYDRAULICS
%
% Whole-plant hydraulic solver with:
%
%   nonlinear soil -> root transport
%   root hydraulic source selection
%   vulnerability-dependent stem and leaf conductance
%   stomatal regulation
%   xylem hydraulic capacitance
%
%
% Hydraulic topology:
%
% soil -> root -> stem -> [xylem storage node] -> leaf -> atmosphere
%
%
% With capacitance:
%
%   Jroot ~= E during transient conditions
%
% and:
%
%   Jroot - E = Cx * dPsiX/dt
%
%
% Positive:
%
%   E - Jroot > 0
%
% means xylem storage is releasing water.
%
%
% REQUIRED INPUT
%
% timeDays
%
% xylemTraits may contain:
%
%   capacitanceEnabled
%
% and either:
%
%   capacitance_mol_m2_MPa
%
% or:
%
%   capacitanceTauDays
%
%
% If capacitance is disabled, the equations reduce to the previous
% steady-state hydraulic formulation.

if nargin < 7 || isempty(leafTraits)

    leafTraits = ...
        struct;

end

if nargin < 8 || isempty(initialState)

    initialState = ...
        struct;

end

if ~isfield(leafTraits,'allowReverseXylemApoplastFlow')
    leafTraits.allowReverseXylemApoplastFlow = false;
end

if ~isfield(leafTraits,'allowReverseRootFlow')
    leafTraits.allowReverseRootFlow = false;
end

if ~isfield(leafTraits,'maxReverseRootFlowFraction')
    leafTraits.maxReverseRootFlowFraction = Inf;
end

%% ------------------------------------------------------------------------
% Input dimensions
% -------------------------------------------------------------------------

timeDays = timeDays(:).';

psiSoil = soil.psiMPa;

[nLayer,nTime] = size(psiSoil);

%% ------------------------------------------------------------------------
% Optional continuation state
%
% This supplies the hydraulic state immediately preceding timeDays(1).
% Existing callers remain unchanged when initialState is empty.
% -------------------------------------------------------------------------

hasInitialState = ...
    isfield(initialState,'timeDays') ...
    && isscalar(initialState.timeDays) ...
    && isfinite(initialState.timeDays);


if hasInitialState

    if initialState.timeDays >= timeDays(1)

        error( ...
            'initialState.timeDays must precede timeDays(1).');

    end

end

if numel(timeDays) ~= nTime
    error('timeDays must match the number of soil timesteps.');
end

if any(diff(timeDays) <= 0)
    error('timeDays must be strictly increasing.');
end


%% ------------------------------------------------------------------------
% Required plant hydraulic parameters
% -------------------------------------------------------------------------

required = { ...
    'KstemMax', ...
    'KleafMax', ...
    'P50Stem', ...
    'P50Leaf', ...
    'vulnerabilityShapeStem', ...
    'vulnerabilityShapeLeaf', ...
    'height_m', ...
    'gsMax', ...
    'psi50Stomata', ...
    'stomatalSlope', ...
    'psiMinMPa', ...
    'psiMaxMPa'};

for i = 1:numel(required)

    if ~isfield(plantState,required{i})

        error( ...
            'plantState.%s is required.', ...
            required{i});

    end

end

if isfield(plantState,'leafParallelHydraulics')

    leafNetCfg = ...
        plantState.leafParallelHydraulics;

else

    leafNetCfg = ...
        struct;


    leafNetCfg.enabled = ...
        false;

end



if leafNetCfg.enabled

    %% ------------------------------------------------------------
    % Xylem -> apoplast conductance mode
    %
    % "fixed"
    %     use leafNetCfg.KxylemApoplast
    %
    % "legacyvulnerability"
    %     use the existing vulnerability-dependent Kleaf relation.
    %     This is required for exact regression against the previous
    %     lumped-leaf solver.
    % ------------------------------------------------------------

    if isfield(leafNetCfg,'KxylemApoplastMode')

        leafNetCfg.KxylemApoplastMode = ...
            lower(string( ...
                leafNetCfg.KxylemApoplastMode));

    else

        leafNetCfg.KxylemApoplastMode = ...
            "fixed";

    end


    validKXAModes = [ ...
        "fixed", ...
        "legacyvulnerability"];


    if ~any( ...
            leafNetCfg.KxylemApoplastMode ...
            == validKXAModes)

        error( ...
            'Unknown leafParallelHydraulics.KxylemApoplastMode.');

    end


    %% ------------------------------------------------------------
    % Conductances always required by parallel network
    % ------------------------------------------------------------

    requiredLeafNet = { ...
        'KapoplastSymplast', ...
        'KsymplastEvaporation', ...
        'KapoplastEvaporation'};


    for q = 1:numel(requiredLeafNet)

        if ~isfield(leafNetCfg,requiredLeafNet{q})

            error( ...
                'plantState.leafParallelHydraulics.%s is required.', ...
                requiredLeafNet{q});

        end

    end


    if leafNetCfg.KapoplastSymplast <= 0

        error( ...
            'KapoplastSymplast must be positive.');

    end


    if leafNetCfg.KsymplastEvaporation <= 0

        error( ...
            'KsymplastEvaporation must be positive.');

    end


    if leafNetCfg.KapoplastEvaporation < 0

        error( ...
            'KapoplastEvaporation must be nonnegative.');

    end


    %% Fixed KXA requires an explicit value.

    if leafNetCfg.KxylemApoplastMode == "fixed"

        if ~isfield(leafNetCfg,'KxylemApoplast')

            error( ...
                'KxylemApoplast is required for fixed KXA mode.');

        end


        if ...
                ~isscalar(leafNetCfg.KxylemApoplast) ...
                || ~isfinite(leafNetCfg.KxylemApoplast) ...
                || leafNetCfg.KxylemApoplast <= 0

            error( ...
                'KxylemApoplast must be positive.');

        end

    end

end

%% ------------------------------------------------------------------------
% Constants
% -------------------------------------------------------------------------

secondsPerDay = 86400;

gravityMPa = ...
    0.00981 ...
    .* plantState.height_m;


%% ========================================================================
% Xylem storage model
% ========================================================================

if isempty(xylemTraits)

    xylemTraits = ...
        struct;

end


if ~isfield(xylemTraits,'storageModel')

    storageModel = ...
        "constantcapacitance";

else

    storageModel = ...
        lower(string( ...
            xylemTraits.storageModel));

end


validStorageModels = [ ...
    "constantcapacitance", ...
    "twocomponentpiecewise", ...
    "twocomponenthysteretic"];


if ~any(storageModel == validStorageModels)

    error( ...
        'Unrecognized xylemTraits.storageModel.');

end


%% ========================================================================
% Xylem capacitance activation
% ========================================================================

if isfield(xylemTraits,'capacitanceEnabled')

    if ~isscalar(xylemTraits.capacitanceEnabled)

        error( ...
            'xylemTraits.capacitanceEnabled must be scalar.');

    end


    capacitanceEnabled = ...
        logical(xylemTraits.capacitanceEnabled);

else

    switch storageModel

        case "constantcapacitance"

            capacitanceEnabled = ...
                isfield( ...
                    xylemTraits, ...
                    'capacitance_mol_m2_MPa') ...
                || ...
                isfield( ...
                    xylemTraits, ...
                    'capacitanceTauDays');


        case { ...
                "twocomponentpiecewise", ...
                "twocomponenthysteretic"}

            % Selecting a pressure-dependent storage model implies
            % dynamic xylem storage unless explicitly disabled.

            capacitanceEnabled = ...
                true;

    end

end


%% ========================================================================
% Diagnostic dynamic-storage-flux ablation
% ========================================================================

diagnosticDisableXylemStorageFlux = ...
    false;

if isfield(xylemTraits,'diagnosticDisableDynamicStorageFlux')

    if ~isscalar(xylemTraits.diagnosticDisableDynamicStorageFlux)
        error([ ...
            'xylemTraits.diagnosticDisableDynamicStorageFlux must be ' ...
            'scalar.']);
    end

    diagnosticDisableXylemStorageFlux = ...
        logical(xylemTraits.diagnosticDisableDynamicStorageFlux);
end

%% ========================================================================
% Xylem capacitance parameterization
% ========================================================================

Cx = ...
    0;


capacitanceParameterization = ...
    'disabled';


if capacitanceEnabled

    switch storageModel


        %% ----------------------------------------------------------------
        % Constant capacitance
        %
        % Requires an explicit constant Cx or a reference time constant.
        % -----------------------------------------------------------------

        case "constantcapacitance"

            if isfield( ...
                    xylemTraits, ...
                    'capacitance_mol_m2_MPa')

                Cx = ...
                    xylemTraits.capacitance_mol_m2_MPa;


                capacitanceParameterization = ...
                    'direct capacitance';


            elseif isfield( ...
                    xylemTraits, ...
                    'capacitanceTauDays')

                Cx = ...
                    plantState.KstemMax ...
                    .* xylemTraits.capacitanceTauDays ...
                    .* secondsPerDay;


                capacitanceParameterization = ...
                    'reference hydraulic time constant';


            else

                error([ ...
                    'Constant-capacitance xylem storage requires ' ...
                    'xylemTraits.capacitance_mol_m2_MPa or ' ...
                    'xylemTraits.capacitanceTauDays.']);

            end


            if ...
                    ~isscalar(Cx) ...
                    || ~isfinite(Cx) ...
                    || Cx <= 0

                error( ...
                    'Xylem capacitance must be positive.');

            end


        %% ----------------------------------------------------------------
        % Pressure-dependent two-component storage
        %
        % Effective capacitance is calculated locally from the
        % pressure-storage relation:
        %
        %   Ceff = dS/dpsiX
        %
        % -----------------------------------------------------------------

        case "twocomponentpiecewise"

            Cx = ...
                NaN;


            capacitanceParameterization = ...
                'two-component pressure-dependent storage';


        case "twocomponenthysteretic"

            Cx = ...
                NaN;


            capacitanceParameterization = ...
                'two-component hysteretic pressure-dependent storage';

    end

end

%% ========================================================================
% Leaf hydraulic storage configuration
%
% Supported constitutive models:
%
%   linear
%   pressureVolume
%
% Leaf storage satisfies:
%
%   dSleaf/dt =
%       JxylemToLeaf - E
%
% The constitutive model determines:
%
%   psiLeaf -> Sleaf
%
% and therefore the effective capacitance:
%
%   Cleaf = dSleaf/dpsiLeaf
% ========================================================================

leafStorageCfg = ...
    struct;


leafStorageCfg.enabled = ...
    false;


if isfield(leafTraits,'storage')

    leafStorageCfg = ...
        leafTraits.storage;


    if ~isfield(leafStorageCfg,'enabled')

        leafStorageCfg.enabled = ...
            false;

    end

end


leafCapacitanceEnabled = ...
    logical(leafStorageCfg.enabled);


% Optional diagnostic-only ablation of the dynamic leaf storage flux.
% The P-V constitutive state remains active because the parallel leaf
% hydraulic network depends on it. Only dSleaf/dt is removed from the
% hydraulic balance. Defaults false.
diagnosticDisableLeafStorageFlux = ...
    false;

if isfield(leafStorageCfg,'diagnosticDisableDynamicStorageFlux')

    if ~isscalar(leafStorageCfg.diagnosticDisableDynamicStorageFlux)
        error([ ...
            'leafTraits.storage.diagnosticDisableDynamicStorageFlux must ' ...
            'be scalar.']);
    end

    diagnosticDisableLeafStorageFlux = ...
        logical(leafStorageCfg.diagnosticDisableDynamicStorageFlux);
end


%% ------------------------------------------------------------------------
% Leaf-storage model
% -------------------------------------------------------------------------

if isfield(leafStorageCfg,'model')

    leafStorageModel = ...
        lower(string(leafStorageCfg.model));

else

    % Preserve historical behavior.
    leafStorageModel = ...
        "linear";

end


validLeafStorageModels = [ ...
    "linear", ...
    "pressurevolume", ...
    "pressure-volume", ...
    "pv"];


if ...
        leafCapacitanceEnabled ...
        && ~any(leafStorageModel == validLeafStorageModels)

    error( ...
        'Unrecognized leafTraits.storage.model.');

end


%% ------------------------------------------------------------------------
% Model-specific validation
% -------------------------------------------------------------------------

if leafCapacitanceEnabled

    switch leafStorageModel


        case "linear"

            requiredLeafStorage = { ...
                'referenceWaterContent_mol_m2', ...
                'referencePsi_MPa', ...
                'capacitance_mol_m2_MPa'};


            for k = 1:numel(requiredLeafStorage)

                if ~isfield( ...
                        leafStorageCfg, ...
                        requiredLeafStorage{k})

                    error( ...
                        'leafTraits.storage.%s is required for linear storage.', ...
                        requiredLeafStorage{k});

                end

            end


            if ...
                    ~isscalar( ...
                        leafStorageCfg.referenceWaterContent_mol_m2) ...
                    || ...
                    ~isfinite( ...
                        leafStorageCfg.referenceWaterContent_mol_m2) ...
                    || ...
                    leafStorageCfg.referenceWaterContent_mol_m2 <= 0

                error( ...
                    'Leaf reference water content must be positive.');

            end


            if ...
                    ~isscalar( ...
                        leafStorageCfg.referencePsi_MPa) ...
                    || ...
                    ~isfinite( ...
                        leafStorageCfg.referencePsi_MPa)

                error( ...
                    'Leaf reference water potential must be finite.');

            end


            if ...
                    ~isscalar( ...
                        leafStorageCfg.capacitance_mol_m2_MPa) ...
                    || ...
                    ~isfinite( ...
                        leafStorageCfg.capacitance_mol_m2_MPa) ...
                    || ...
                    leafStorageCfg.capacitance_mol_m2_MPa < 0

                error( ...
                    'Leaf capacitance must be nonnegative.');

            end


        case { ...
                "pressurevolume", ...
                "pressure-volume", ...
                "pv"}

            requiredLeafStorage = { ...
                'saturatedWaterContent_mol_m2', ...
                'osmoticPotentialFullTurgor_MPa', ...
                'turgorLossPoint_MPa', ...
                'relativeWaterContentTurgorLoss', ...
                'apoplasticFraction'};


            for k = 1:numel(requiredLeafStorage)

                if ~isfield( ...
                        leafStorageCfg, ...
                        requiredLeafStorage{k})

                    error( ...
                        'leafTraits.storage.%s is required for pressure-volume storage.', ...
                        requiredLeafStorage{k});

                end

            end


            % Run the constitutive function once here so its more
            % detailed physiological input checks fail immediately.

            leaf_storage_pressure_volume_state( ...
                leafStorageCfg.turgorLossPoint_MPa, ...
                leafStorageCfg);

    end

end

%% ========================================================================
% Parallel leaf hydraulics currently requires pressure-volume storage
%
% The first two-domain implementation treats:
%
%   apoplast:
%       hydraulically massless / fixed water amount
%
%   symplast:
%       finite P-V storage reservoir
%
% Therefore the parallel solver currently requires the validated
% pressure-volume storage model.
% ========================================================================

if leafNetCfg.enabled

    if ~leafCapacitanceEnabled

        error([ ...
            'Parallel leaf hydraulics currently requires ' ...
            'leaf hydraulic storage to be enabled.']);

    end


    if ~any(leafStorageModel == [ ...
            "pressurevolume", ...
            "pressure-volume", ...
            "pv"])

        error([ ...
            'Parallel leaf hydraulics currently requires the ' ...
            'pressure-volume leaf storage model.']);

    end

end

%% ------------------------------------------------------------------------
% Hydraulic storage reservoir
%
% This is the absolute amount of capacitively accessible stored water.
%
% The current linear capacitance equations determine DeltaS from:
%
%   DeltaS = Cx * DeltaPsiX
%
% hydraulicStorageInitial_mol_m2 then converts this to:
%
%   S_hydraulic(t) = S0 + DeltaS(t)
%
% At this stage S0 is used for bookkeeping and depletion diagnostics.
% -------------------------------------------------------------------------

if isfield(xylemTraits,'hydraulicStorageInitial_mol_m2')

    hydraulicStorageInitial = ...
        xylemTraits.hydraulicStorageInitial_mol_m2;

else

    hydraulicStorageInitial = NaN;

end

%% ------------------------------------------------------------------------
% Fully hydrated / anatomical storage capacity
% -------------------------------------------------------------------------

hydraulicStorageCapacity = ...
    NaN;


switch storageModel


    case "constantcapacitance"

       
        % Preserve the historical configured value if supplied.

        if isfield( ...
                xylemTraits, ...
                'hydraulicStorageCapacity_mol_m2')

            hydraulicStorageCapacity = ...
                xylemTraits.hydraulicStorageCapacity_mol_m2;


        elseif isfinite(hydraulicStorageInitial)

            hydraulicStorageCapacity = ...
                hydraulicStorageInitial;

        end


    case { ...
            "twocomponentpiecewise", ...
            "twocomponenthysteretic"}

        % The two-component constitutive model explicitly defines
        % fully hydrated lumen and living-tissue inventories.

        requiredCapacityFields = { ...
            'storageLumenWaterInitial_mol_m2', ...
            'livingStorageWaterInitial_mol_m2'};


        for k = 1:numel(requiredCapacityFields)

            if ~isfield( ...
                    xylemTraits, ...
                    requiredCapacityFields{k})

                error( ...
                    'xylemTraits.%s is required for two-component storage.', ...
                    requiredCapacityFields{k});

            end

        end


        hydraulicStorageCapacity = ...
            xylemTraits.storageLumenWaterInitial_mol_m2 ...
            + ...
            xylemTraits.livingStorageWaterInitial_mol_m2;

end


if ...
        isfinite(hydraulicStorageCapacity) ...
        && hydraulicStorageCapacity <= 0

    error( ...
        'Xylem hydraulic storage capacity must be positive.');

end

%% ========================================================================
% Xylem hydraulic-storage bookkeeping
% ========================================================================

hydraulicStorageInitial = ...
    NaN;


if isfield( ...
        xylemTraits, ...
        'hydraulicStorageInitial_mol_m2')

    hydraulicStorageInitial = ...
        xylemTraits.hydraulicStorageInitial_mol_m2;

end


%% ========================================================================
% Fully hydrated xylem storage capacity
%
% Constant-capacitance model:
%
%   absolute capacity is not defined by Cx itself. Preserve an explicitly
%   supplied capacity, or historically use the configured initial amount.
%
% Two-component models:
%
%   component maximum inventories define anatomical capacity directly:
%
%       Smax =
%           S_lumen,max
%           + S_living,max
% ========================================================================

hydraulicStorageCapacity = ...
    NaN;


switch storageModel


    case "constantcapacitance"

        if isfield( ...
                xylemTraits, ...
                'hydraulicStorageCapacity_mol_m2')

            hydraulicStorageCapacity = ...
                xylemTraits.hydraulicStorageCapacity_mol_m2;


        elseif isfinite(hydraulicStorageInitial)

            hydraulicStorageCapacity = ...
                hydraulicStorageInitial;

        end


    case { ...
            "twocomponentpiecewise", ...
            "twocomponenthysteretic"}

        requiredCapacityFields = { ...
            'storageLumenWaterInitial_mol_m2', ...
            'livingStorageWaterInitial_mol_m2'};


        for k = 1:numel(requiredCapacityFields)

            if ~isfield( ...
                    xylemTraits, ...
                    requiredCapacityFields{k})

                error( ...
                    'xylemTraits.%s is required for two-component storage.', ...
                    requiredCapacityFields{k});

            end

        end


        hydraulicStorageCapacity = ...
            xylemTraits.storageLumenWaterInitial_mol_m2 ...
            + ...
            xylemTraits.livingStorageWaterInitial_mol_m2;

end


if ...
        isfinite(hydraulicStorageCapacity) ...
        && hydraulicStorageCapacity <= 0

    error( ...
        'Xylem hydraulic storage capacity must be positive.');

end


%% ========================================================================
% Storage-model-specific initial-state validation
% ========================================================================

if capacitanceEnabled

    switch storageModel


        case "constantcapacitance"

            if ...
                    ~isscalar(hydraulicStorageInitial) ...
                    || ~isfinite(hydraulicStorageInitial) ...
                    || hydraulicStorageInitial <= 0

                error([ ...
                    'Constant-capacitance xylem storage requires ' ...
                    'xylemTraits.hydraulicStorageInitial_mol_m2 ' ...
                    'to be positive.']);

            end


        case { ...
                "twocomponentpiecewise", ...
                "twocomponenthysteretic"}

            % No configured initial amount is needed.
            %
            % The actual initial inventory is established from the
            % solved initial xylem potential:
            %
            %   S0 = S(psiX0)

    end

end

%% ------------------------------------------------------------------------
% Prepare nonlinear soil-root transport model
% -------------------------------------------------------------------------

psiMFPMin = min([ ...
    plantState.psiMinMPa, ...
    min(soil.psiMPa,[],'all','omitnan')]);


rootPath = prepare_root_pathway_model( ...
    rootTraits, ...
    soil.hydro, ...
    plantState, ...
    psiMFPMin, ...
    nLayer);


%% ------------------------------------------------------------------------
% Allocate primary states
% -------------------------------------------------------------------------

psiRoot = NaN(1,nTime);
psiLeaf = NaN(1,nTime);
psiXylem = NaN(1,nTime);

%% ------------------------------------------------------------------------
% Two-domain leaf hydraulic states
% -------------------------------------------------------------------------

psiSymplast = ...
    NaN(1,nTime);


psiApoplast = ...
    NaN(1,nTime);


psiEvaporationSite = ...
    NaN(1,nTime);


JxylemToApoplast = ...
    NaN(1,nTime);


JapoplastToSymplast = ...
    NaN(1,nTime);


JsymplastToEvaporation = ...
    NaN(1,nTime);


JapoplastToEvaporation = ...
    NaN(1,nTime);


symplasticAdvectiveFractionHydraulic = ...
    NaN(1,nTime);


apoplasticAdvectiveFractionHydraulic = ...
    NaN(1,nTime);


leafParallelApoplastClosure = ...
    NaN(1,nTime);


leafParallelSymplastClosure = ...
    NaN(1,nTime);


leafParallelEvaporationClosure = ...
    NaN(1,nTime);


leafParallelWholeClosure = ...
    NaN(1,nTime);


leafParallelStorageRateResidual = ...
    NaN(1,nTime);


KxylemApoplastOut = ...
    NaN(1,nTime);

E  = NaN(1,nTime);
gs = NaN(1,nTime);

rootFlux = NaN(nLayer,nTime);
rootWeights = NaN(nLayer,nTime);

KstemOut = NaN(1,nTime);
KleafOut = NaN(1,nTime);

rootUptakeFlux = ...
    NaN(nLayer,nTime);


rootReleaseFlux = ...
    NaN(nLayer,nTime);


totalRootUptakeFlux = ...
    NaN(1,nTime);


totalRootReleaseFlux = ...
    NaN(1,nTime);

%% ------------------------------------------------------------------------
% Capacitance diagnostics
% -------------------------------------------------------------------------

requiredStemFlow = NaN(1,nTime);
actualStemFlow = NaN(1,nTime);

storageRate = NaN(1,nTime);

netStorageReleaseFlux = NaN(1,nTime);
storageReleaseFlux = NaN(1,nTime);
storageRechargeFlux = NaN(1,nTime);

storageContributionFraction = NaN(1,nTime);

deltaXylemStorage = NaN(1,nTime);
xylemStorageWaterDeficit = NaN(1,nTime);

capacitanceMassBalanceResidual = NaN(1,nTime);
stemFlowResidual = NaN(1,nTime);

%% ------------------------------------------------------------------------
% Pressure-dependent storage diagnostics
% -------------------------------------------------------------------------

hydraulicStorageAmount = ...
    NaN(1,nTime);

hydraulicStorageLumenAmount = ...
    NaN(1,nTime);

hydraulicStorageLivingAmount = ...
    NaN(1,nTime);

xylemCapacitanceEffective = ...
    NaN(1,nTime);

xylemCapacitanceLumen = ...
    NaN(1,nTime);

xylemCapacitanceLiving = ...
    NaN(1,nTime);

%% ------------------------------------------------------------------------
% General solver diagnostics
% -------------------------------------------------------------------------

hydraulicResidualOut = NaN(1,nTime);

solverFallback = false(1,nTime);

nHydraulicRoots = NaN(1,nTime);

allHydraulicRoots = cell(1,nTime);

% Detailed diagnostics are populated only at timesteps where the outer
% hydraulic root cannot be located. Successful timesteps retain [].
failureDiagnostics = cell(1,nTime);

branchGs = cell(1,nTime);
branchE = cell(1,nTime);
branchKleaf = cell(1,nTime);

hydraulicallyFeasible = false(1,nTime);

massBalanceResidual = NaN(1,nTime);


%% ------------------------------------------------------------------------
% Soil-root diagnostics
% -------------------------------------------------------------------------

KrootTissueOut = NaN(nLayer,nTime);

KsoilPathOut = NaN(nLayer,nTime);

KeffectiveOut = NaN(nLayer,nTime);

soilLimitOut = NaN(nLayer,nTime);

rootLimitOut = NaN(nLayer,nTime);

soilKrelOut = NaN(nLayer,nTime);

KsoilUnsatOut = NaN(nLayer,nTime);

KsoilPathBulkOut = NaN(nLayer,nTime);

soilFluxResidualOut = NaN(nLayer,nTime);

psiRootSurfaceOut = NaN(nLayer,nTime);

storageLumenBranchCode = NaN(1,nTime);

%% ------------------------------------------------------------------------
% Leaf hydraulic storage
% -------------------------------------------------------------------------

leafHydraulicStorageAmount = ...
    NaN(1,nTime);


leafHydraulicStorageRate = ...
    zeros(1,nTime);


JxylemToLeaf = ...
    NaN(1,nTime);


leafStorageReleaseFlux = ...
    zeros(1,nTime);


leafStorageRechargeFlux = ...
    zeros(1,nTime);


leafStorageNetReleaseFlux = ...
    zeros(1,nTime);


leafStorageFluxClosureResidual = ...
    NaN(1,nTime);


leafCapacitanceEffective = ...
    NaN(1,nTime);


leafStorageBoundActive = ...
    false(1,nTime);


leafHydraulicStorageInitialActual = ...
    NaN;

leafRWC = ...
    NaN(1,nTime);


leafSymplasticRWC = ...
    NaN(1,nTime);


leafOsmoticPsi = ...
    NaN(1,nTime);


leafTurgorPsi = ...
    NaN(1,nTime);


leafApoplasticWater = ...
    NaN(1,nTime);


leafSymplasticWater = ...
    NaN(1,nTime);

%% ------------------------------------------------------------------------
% Component-specific hydraulic storage fluxes
% -------------------------------------------------------------------------

storageLumenRate = ...
    NaN(1,nTime);

storageLivingRate = ...
    NaN(1,nTime);


storageLumenReleaseFlux = ...
    NaN(1,nTime);

storageLumenRechargeFlux = ...
    NaN(1,nTime);


storageLivingReleaseFlux = ...
    NaN(1,nTime);

storageLivingRechargeFlux = ...
    NaN(1,nTime);


componentStorageRateResidual = ...
    NaN(1,nTime);

componentStorageNetFluxResidual = ...
    NaN(1,nTime);

usedContinuation = ...
    false(1,nTime);


continuationResidual = ...
    NaN(1,nTime);


continuationStateValid = ...
    false(1,nTime);

%% ------------------------------------------------------------------------
% Solve each timestep
% -------------------------------------------------------------------------

psiX0 = NaN;
hydraulicStorageInitialActual = NaN;
residualEvalCount = 0;

for t = 1:nTime
    
    residualEvalCount = residualEvalCount + 1;

    psiS = ...
        psiSoil(:,t);


    %% ------------------------------------------------------------
    % Dynamic state entering timestep
    % ------------------------------------------------------------

    if t > 1

        dtSeconds = ...
            ( ...
            timeDays(t) ...
            - timeDays(t-1) ...
            ) ...
            .* secondsPerDay;

    elseif hasInitialState

        dtSeconds = ...
            ( ...
            timeDays(1) ...
            - initialState.timeDays ...
            ) ...
            .* secondsPerDay;

    else

        dtSeconds = ...
            NaN;

    end


    %% Xylem-storage state

    if ...
            capacitanceEnabled ...
            && t > 1 ...
            && isfinite(psiXylem(t-1))

        psiXPrevious = ...
            psiXylem(t-1);

        useCapacitance = ...
            ~diagnosticDisableXylemStorageFlux;

    elseif ...
            capacitanceEnabled ...
            && t == 1 ...
            && hasInitialState ...
            && isfield(initialState,'psiXylem') ...
            && isfinite(initialState.psiXylem)

        psiXPrevious = ...
            initialState.psiXylem;

        useCapacitance = ...
            ~diagnosticDisableXylemStorageFlux;

    else

        useCapacitance = ...
            false;

        psiXPrevious = ...
            NaN;

    end


    %% Leaf-storage state

    if ...
            leafCapacitanceEnabled ...
            && t > 1 ...
            && isfinite(leafHydraulicStorageAmount(t-1))

        leafStoragePrevious = ...
            leafHydraulicStorageAmount(t-1);

        useLeafCapacitance = ...
            ~diagnosticDisableLeafStorageFlux;

    elseif ...
            leafCapacitanceEnabled ...
            && t == 1 ...
            && hasInitialState ...
            && isfield(initialState,'leafHydraulicStorageAmount_mol_m2') ...
            && isfinite(initialState.leafHydraulicStorageAmount_mol_m2)

        leafStoragePrevious = ...
            initialState.leafHydraulicStorageAmount_mol_m2;

        useLeafCapacitance = ...
            ~diagnosticDisableLeafStorageFlux;

    else

        useLeafCapacitance = ...
            false;

        leafStoragePrevious = ...
            NaN;

    end

    %% Previous symplastic water state

    if ...
            leafNetCfg.enabled ...
            && t > 1 ...
            && isfinite(leafSymplasticWater(t-1))

        leafSymplasticWaterPrevious = ...
            leafSymplasticWater(t-1);

    elseif ...
            leafNetCfg.enabled ...
            && t == 1 ...
            && hasInitialState ...
            && isfield(initialState,'leafSymplasticWaterAmount_mol_m2') ...
            && isfinite(initialState.leafSymplasticWaterAmount_mol_m2)

        leafSymplasticWaterPrevious = ...
            initialState.leafSymplasticWaterAmount_mol_m2;

    else

        leafSymplasticWaterPrevious = ...
            NaN;

    end


    %% ------------------------------------------------------------
    % Coupled hydraulic residual
    % ------------------------------------------------------------
    if ...
            t > 1 ...
            && isfinite(hydraulicStorageAmount(t-1))

        storagePrevious = ...
            hydraulicStorageAmount(t-1);

    elseif ...
            t == 1 ...
            && hasInitialState ...
            && isfield(initialState,'hydraulicStorageAmount_mol_m2') ...
            && isfinite(initialState.hydraulicStorageAmount_mol_m2)

        storagePrevious = ...
            initialState.hydraulicStorageAmount_mol_m2;

    else

        storagePrevious = ...
            NaN;

    end
    if ...
            t > 1 ...
            && isfinite(hydraulicStorageLumenAmount(t-1))

        storageLumenPrevious = ...
            hydraulicStorageLumenAmount(t-1);

    elseif ...
            t == 1 ...
            && hasInitialState ...
            && isfield(initialState,'hydraulicStorageLumenAmount_mol_m2') ...
            && isfinite(initialState.hydraulicStorageLumenAmount_mol_m2)

        storageLumenPrevious = ...
            initialState.hydraulicStorageLumenAmount_mol_m2;

    else

        storageLumenPrevious = ...
            NaN;

    end

    % Search-only candidate traits.
    % Reverse-flow states are permitted only while evaluating/bracketing
    % the nonlinear residual. The final selected hydraulic state is
    % re-evaluated below using the original strict leafTraits.

    leafTraitsSearch = ...
        leafTraits;

    leafTraitsSearch.allowReverseXylemApoplastFlow = ...
        true;

    leafTraitsSearch.allowReverseRootFlow = ...
        true;

    leafTraitsSearch.maxReverseRootFlowFraction = ...
        Inf;

    % Root-supply branch continuation reference.
    %
    % Under strong drought the nonlinear soil-root supply curve can contain
    % multiple admissible roots for the same Jrequired. Follow the branch
    % nearest the previously accepted root pressure instead of imposing an
    % arbitrary least-negative-root rule.
    if ...
            t > 1 ...
            && isfinite(psiRoot(t-1))

        psiRootBranchReference = ...
            psiRoot(t-1);

    elseif ...
            t == 1 ...
            && hasInitialState ...
            && isfield(initialState,'psiRoot') ...
            && isfinite(initialState.psiRoot)

        psiRootBranchReference = ...
            initialState.psiRoot;

    else

        psiRootBranchReference = ...
            NaN;

    end

    % Precompute the nonlinear soil-root supply curve ONCE for this soil
    % state. The curve depends on psiSoil/rootPath, not on psiLeaf.
    %
    if isfield(rootPath,'rootDemandPsiMinMPa')
        rootCachePsiMin = rootPath.rootDemandPsiMinMPa;
    else
        rootCachePsiMin = -5;
    end

    rootDemandBracketCache = ...
        buildRootDemandBracketCacheLocal( ...
            psiS, ...
            rootPath, ...
            rootCachePsiMin, ...
            0, ...
            33);

    fun = @(psiL) ...
    hydraulicResidualCapacitance( ...
        psiL, ...
        psiS, ...
        rootPath, ...
        met, ...
        plantState, ...
        t, ...
        psiXPrevious, ...
        storagePrevious, ...
        storageLumenPrevious, ...
        dtSeconds, ...
        Cx, ...
        useCapacitance, ...
        storageModel, ...
        xylemTraits, ...
        leafStoragePrevious, ...
        useLeafCapacitance, ...
        leafStorageCfg, ...
        leafNetCfg, ...
        leafSymplasticWaterPrevious, ...
        gravityMPa, ...
        leafTraitsSearch, ...
        rootDemandBracketCache, ...
        psiRootBranchReference);


    %% ------------------------------------------------------------
    % Previous-state continuation guess
    % ------------------------------------------------------------

    if ...
            t > 1 ...
            && isfinite(psiLeaf(t-1))

        psiGuess = ...
            psiLeaf(t-1);

    elseif ...
            t == 1 ...
            && hasInitialState ...
            && isfield(initialState,'psiLeaf') ...
            && isfinite(initialState.psiLeaf)

        psiGuess = ...
            initialState.psiLeaf;

    else

        psiGuess = ...
            -1.0;

    end


    psiGuess = ...
        min( ...
        max( ...
        psiGuess, ...
        plantState.psiMinMPa), ...
        plantState.psiMaxMPa);


    %% ------------------------------------------------------------
    % Previous-state continuation
    %
    % Dynamic storage can make the admissible hydraulic-root region
    % extremely narrow at short timesteps. Test the previous accepted
    % leaf state explicitly before performing the global root search.
    % ------------------------------------------------------------

    [residualGuess,stateGuess] = ...
        fun(psiGuess);


    continuationResidual(t) = ...
        residualGuess;


    continuationStateValid(t) = ...
        stateGuess.valid ...
        && isfinite(residualGuess);


    continuationTolerance = ...
        1e-10;


    %% ------------------------------------------------------------
    % Locate physical hydraulic root
    % ------------------------------------------------------------

    if ...
            (t > 1 || hasInitialState) ...
            && continuationStateValid(t) ...
            && abs(residualGuess) <= continuationTolerance

        % Previous accepted state still satisfies the hydraulic equations.

        psiL = ...
            psiGuess;


        solverInfo = ...
            struct;


        solverInfo.rootFound = ...
            true;


        solverInfo.nRoots = ...
            1;


        solverInfo.allRoots = ...
            psiL;


        solverInfo.usedContinuation = ...
            true;

    else

        % Previous state is not itself an exact root.
        % Fall back to the existing global search.

        [psiL,solverInfo] = ...
            findPhysicalHydraulicRoot( ...
            fun, ...
            plantState.psiMinMPa, ...
            plantState.psiMaxMPa, ...
            psiGuess);


        solverInfo.usedContinuation = ...
            false;

    end


    usedContinuation(t) = ...
        solverInfo.usedContinuation;


    solverFallback(t) = ...
        ~solverInfo.rootFound;


    nHydraulicRoots(t) = ...
        solverInfo.nRoots;


    allHydraulicRoots{t} = ...
        solverInfo.allRoots;


    %% ------------------------------------------------------------
    % Failure diagnostics
    % ------------------------------------------------------------

    if ~solverInfo.rootFound

        fprintf('\nHYDRAULIC FAILURE AT INDEX = %d\n',t);

        fprintf( ...
            'Time = %.3f days\n', ...
            timeDays(t));


        psiDebug = linspace( ...
            plantState.psiMinMPa, ...
            plantState.psiMaxMPa, ...
            5001);


        Fdebug = NaN(size(psiDebug));
        rejectionCodeDebug = strings(size(psiDebug));
        for q = 1:numel(psiDebug)

            [Fdebug(q),state_q] = ...
                fun(psiDebug(q));

            if isfield(state_q,'rejectCode')
                rejectionCodeDebug(q) = string(state_q.rejectCode);
            else
                rejectionCodeDebug(q) = "unclassified";
            end

        end


        finiteDebug = ...
            isfinite(Fdebug);


        fprintf( ...
            'Finite residual points: %d / %d\n', ...
            sum(finiteDebug), ...
            numel(Fdebug));


        closestIndex = NaN;
        closestState = struct;


        if any(finiteDebug)

            fprintf( ...
                'Residual min = %.8g\n', ...
                min(Fdebug(finiteDebug)));

            fprintf( ...
                'Residual max = %.8g\n', ...
                max(Fdebug(finiteDebug)));


            [~,qmin] = ...
                min(abs(Fdebug(finiteDebug)));

            ii = ...
                find(finiteDebug);

            closestIndex = ii(qmin);
            [~,closestState] = fun(psiDebug(closestIndex));


            fprintf( ...
                'Closest residual to zero = %.8g\n', ...
                Fdebug(closestIndex));

            fprintf( ...
                'At psiLeaf = %.8f MPa\n', ...
                psiDebug(closestIndex));

            if closestState.valid

                fprintf('\n--- FAILURE STATE ---\n')

                fprintf('psiRoot      = %.6e MPa\n',closestState.psiRoot)
                fprintf('psiXylem     = %.6e MPa\n',closestState.psiXylem)
                fprintf('psiLeaf      = %.6e MPa\n',psiDebug(closestIndex))

                fprintf('Jroot        = %.6e\n',closestState.Jroot)
                fprintf('Jrequired    = %.6e\n',closestState.requiredStemFlow)
                fprintf('Jstem        = %.6e\n',closestState.actualStemFlow)

                fprintf('root residual = %.6e\n', ...
                    closestState.Jroot - closestState.requiredStemFlow)

                fprintf('stem margin = %.6e\n', ...
                    closestState.actualStemFlow - closestState.requiredStemFlow)

            end

        end


        [reasonNames,~,reasonIndex] = unique(rejectionCodeDebug);
        reasonCounts = accumarray(reasonIndex(:),1);

        [reasonCounts,order] = sort(reasonCounts,'descend');
        reasonNames = reasonNames(order);

        fprintf('Rejection summary:\n');

        for q = 1:numel(reasonNames)
            fprintf('  %-34s %4d\n',char(reasonNames(q)),reasonCounts(q));
        end

        %% ------------------------------------------------------------
        % Enumerate ALL root-demand branches at this failed timestep
        % ------------------------------------------------------------
        %
        % The soil-root supply curve depends on psiS and rootPath but not
        % on psiLeaf. Evaluate it once, then intersect it with each valid
        % leaf candidate's Jrequired. This reveals whether the outer stem
        % residual could cross zero on a root branch other than the one
        % selected by the production inner solve.

        fprintf('\nROOT-BRANCH ENUMERATION AT FAILED TIMESTEP\n');

        psiRootBranchGrid = ...
            linspace(-5,0,101);

        rootSupplyGrid = ...
            NaN(size(psiRootBranchGrid));

        for rri = 1:numel(psiRootBranchGrid)

            LL = ...
                soil_root_layer_flux( ...
                    psiS, ...
                    psiRootBranchGrid(rri), ...
                    rootPath);

            rootSupplyGrid(rri) = ...
                sum(LL.flux,'omitnan');

        end

        branchPsi = ...
            NaN(numel(psiDebug),6);

        branchStemResidual = ...
            NaN(numel(psiDebug),6);

        branchJrequired = ...
            NaN(numel(psiDebug),1);

        branchPsiX = ...
            NaN(numel(psiDebug),1);

        nRootBranches = ...
            zeros(numel(psiDebug),1);

        for q = 1:numel(psiDebug)

            if ~finiteDebug(q)
                continue
            end

            [~,sq] = ...
                fun(psiDebug(q));

            if ...
                    ~isfield(sq,'requiredStemFlow') ...
                    || ~isfinite(sq.requiredStemFlow) ...
                    || ~isfield(sq,'psiXylem') ...
                    || ~isfinite(sq.psiXylem)
                continue
            end

            Jreq_q = ...
                sq.requiredStemFlow;

            psiX_q = ...
                sq.psiXylem;

            rr = ...
                rootSupplyGrid - Jreq_q;

            roots_q = [];

            % Exact sampled roots.
            iex = ...
                find(isfinite(rr) & abs(rr)<=1e-10);

            if ~isempty(iex)
                roots_q = [roots_q psiRootBranchGrid(iex)]; 
            end

            % Sign-changing intervals.
            for kk = 1:(numel(rr)-1)

                if ...
                        isfinite(rr(kk)) ...
                        && isfinite(rr(kk+1)) ...
                        && sign(rr(kk)) ~= sign(rr(kk+1))

                    b = ...
                        sort([ ...
                            psiRootBranchGrid(kk), ...
                            psiRootBranchGrid(kk+1)]);

                    rootFun_q = @(pr) ...
                        sum(soil_root_layer_flux( ...
                            psiS, ...
                            pr, ...
                            rootPath).flux,'omitnan') ...
                        - Jreq_q;

                    try
                        rq = fzero(rootFun_q,b);
                        roots_q(end+1) = rq; 
                    catch
                    end

                end

            end

            if isempty(roots_q)
                continue
            end

            % Deduplicate and order least-negative -> most-negative.
            roots_q = ...
                unique(round(roots_q,10));

            roots_q = ...
                sort(roots_q,'descend');

            nr = ...
                min(numel(roots_q),size(branchPsi,2));

            nRootBranches(q) = ...
                numel(roots_q);

            branchJrequired(q) = ...
                Jreq_q;

            branchPsiX(q) = ...
                psiX_q;

            for kk = 1:nr

                pr = ...
                    roots_q(kk);

                Kstem_q = ...
                    hydraulic_vulnerability( ...
                        pr, ...
                        plantState.KstemMax, ...
                        plantState.P50Stem, ...
                        plantState.vulnerabilityShapeStem);

                Jstem_q = ...
                    Kstem_q ...
                    .* (pr - psiX_q - gravityMPa);

                branchPsi(q,kk) = ...
                    pr;

                branchStemResidual(q,kk) = ...
                    Jstem_q - Jreq_q;

            end

        end

        fprintf( ...
            'Valid leaf candidates enumerated: %d\n', ...
            sum(any(isfinite(branchPsi),2)));

        fprintf( ...
            'Candidates with >1 root branch:    %d\n', ...
            sum(nRootBranches>1));

        fprintf( ...
            'Maximum root branches:   %d\n', ...
            max(nRootBranches));

        allBranchResidual = ...
            branchStemResidual(isfinite(branchStemResidual));

        if ~isempty(allBranchResidual)

            fprintf( ...
                'All-branch stem residual min:   %.10g\n', ...
                min(allBranchResidual));

            fprintf( ...
                'All-branch stem residual max:   %.10g\n', ...
                max(allBranchResidual));

            [~,linClosest] = ...
                min(abs(branchStemResidual(:)));

            [qClosest,bClosest] = ...
                ind2sub(size(branchStemResidual),linClosest);

            fprintf( ...
                'Closest all-branch residual:    %.10g\n', ...
                branchStemResidual(qClosest,bClosest));

            fprintf( ...
                '  psiLeaf:                      %.10f MPa\n', ...
                psiDebug(qClosest));

            fprintf( ...
                '  psiRoot branch:               %.10f MPa\n', ...
                branchPsi(qClosest,bClosest));

            fprintf( ...
                '  psiXylem:                     %.10f MPa\n', ...
                branchPsiX(qClosest));

            fprintf( ...
                '  Jrequired:                    %.10g\n', ...
                branchJrequired(qClosest));

        end

        % Check each root-rank branch for an outer stem-residual sign change
        % across adjacent leaf-pressure candidates.
        anyBranchCrossing = false;

        for bb = 1:size(branchStemResidual,2)

            yb = ...
                branchStemResidual(:,bb);

            ib = ...
                find(isfinite(yb));

            if numel(ib)<2
                continue
            end

            for kk = 1:(numel(ib)-1)

                i1 = ib(kk);
                i2 = ib(kk+1);

                if ...
                        i2 == i1+1 ...
                        && sign(yb(i1)) ~= sign(yb(i2))

                    fprintf( ...
                        ['BRANCH %d OUTER SIGN CHANGE: psiLeaf %.10f -> ' ...
                         '%.10f MPa; residual %.10g -> %.10g\n'], ...
                        bb, ...
                        psiDebug(i1), ...
                        psiDebug(i2), ...
                        yb(i1), ...
                        yb(i2));

                    anyBranchCrossing = true;

                end

            end

        end

        fprintf( ...
            'Any resolved branch outer crossing: %d\n', ...
            anyBranchCrossing);

        fprintf('END ROOT-BRANCH ENUMERATION\n\n');


        FD = struct;
        FD.timeIndex = t;
        FD.timeDays = timeDays(t);
        FD.psiLeafGrid_MPa = psiDebug;
        FD.residual = Fdebug;
        FD.rejectionCode = rejectionCodeDebug;
        FD.reasonNames = reasonNames;
        FD.reasonCounts = reasonCounts;
        FD.nFiniteResidual = sum(finiteDebug);
        FD.closestIndex = closestIndex;

        if isfinite(closestIndex)
            FD.closestPsiLeaf_MPa = psiDebug(closestIndex);
            FD.closestResidual = Fdebug(closestIndex);
            FD.closestState = closestState;
        else
            FD.closestPsiLeaf_MPa = NaN;
            FD.closestResidual = NaN;
            FD.closestState = struct;
        end

        failureDiagnostics{t} = FD;


        % Stop cleanly at first failed timestep

        break

    end


    %% ------------------------------------------------------------
    % Branch diagnostics
    % ------------------------------------------------------------

    roots_t = ...
        solverInfo.allRoots;

    gs_t_branch = NaN(size(roots_t));
    E_t_branch  = NaN(size(roots_t));
    K_t_branch  = NaN(size(roots_t));


    for j = 1:numel(roots_t)

        psiBranch = ...
            roots_t(j);


        [Ebranch,gsBranch] = ...
            transpirationDemand( ...
                psiBranch, ...
                met, ...
                plantState, ...
                t);


        if leafNetCfg.enabled

            Kbranch = ...
                resolveLeafXylemApoplastConductance( ...
                psiBranch, ...
                plantState, ...
                leafNetCfg);

        else

            Kbranch = ...
                hydraulic_vulnerability( ...
                psiBranch, ...
                plantState.KleafMax, ...
                plantState.P50Leaf, ...
                plantState.vulnerabilityShapeLeaf);

        end


        gs_t_branch(j) = gsBranch;
        E_t_branch(j)  = Ebranch;
        K_t_branch(j)  = Kbranch;

    end


    branchGs{t} = ...
        gs_t_branch;

    branchE{t} = ...
        E_t_branch;

    branchKleaf{t} = ...
        K_t_branch;


    %% ------------------------------------------------------------
    % Evaluate selected hydraulic state
    % ------------------------------------------------------------

    [residual_t,state] = ...
    hydraulicResidualCapacitance( ...
        psiL, ...
        psiS, ...
        rootPath, ...
        met, ...
        plantState, ...
        t, ...
        psiXPrevious, ...
        storagePrevious, ...
        storageLumenPrevious, ...
        dtSeconds, ...
        Cx, ...
        useCapacitance, ...
        storageModel, ...
        xylemTraits, ...
        leafStoragePrevious, ...
        useLeafCapacitance, ...
        leafStorageCfg, ...
        leafNetCfg, ...
        leafSymplasticWaterPrevious, ...
        gravityMPa, ...
        leafTraits, ...
        rootDemandBracketCache, ...
        psiRootBranchReference);


    hydraulicResidualOut(t) = ...
        residual_t;


    if ...
            ~state.valid ...
            || ~isfinite(state.psiRoot) ...
            || isempty(state.layer)

        solverFallback(t) = true;

        break

    end


    %% ------------------------------------------------------------
    % Extract state
    % ------------------------------------------------------------

    E_t = ...
        state.E;

    gs_t = ...
        state.gs;

    psiR = ...
        state.psiRoot;

    psiX = ...
        state.psiXylem;

    layer = ...
        state.layer;

    Kstem_t = ...
        state.Kstem;

    Kleaf_t = ...
        state.Kleaf;

    JxylemToLeaf_t = ...
        state.JxylemToLeaf;

    if leafNetCfg.enabled

        psiSymplast(t) = ...
            state.psiSymplast;


        psiApoplast(t) = ...
            state.psiApoplast;


        psiEvaporationSite(t) = ...
            state.psiEvaporationSite;


        JxylemToApoplast(t) = ...
            state.JxylemToApoplast;


        JapoplastToSymplast(t) = ...
            state.JapoplastToSymplast;


        JsymplastToEvaporation(t) = ...
            state.JsymplastToEvaporation;


        JapoplastToEvaporation(t) = ...
            state.JapoplastToEvaporation;


        symplasticAdvectiveFractionHydraulic(t) = ...
            state.symplasticAdvectiveFraction;


        apoplasticAdvectiveFractionHydraulic(t) = ...
            state.apoplasticAdvectiveFraction;


        leafParallelApoplastClosure(t) = ...
            state.leafParallelApoplastClosure;


        leafParallelSymplastClosure(t) = ...
            state.leafParallelSymplastClosure;


        leafParallelEvaporationClosure(t) = ...
            state.leafParallelEvaporationClosure;


        leafParallelWholeClosure(t) = ...
            state.leafParallelWholeClosure;


        leafParallelStorageRateResidual(t) = ...
            state.leafParallelStorageRateResidual;


        KxylemApoplastOut(t) = ...
            state.KxylemApoplast;

    else

        % Legacy compatibility.

        psiSymplast(t) = ...
            psiL;

    end

    leafStorageRate_t = ...
        state.leafStorageRate;


    leafStorageAmount_t = ...
        state.leafStorageAmount;


    leafCapacitanceEffective_t = ...
        state.leafCapacitanceEffective;


    leafStorageBound_t = ...
        state.leafStorageBoundActive;

    leafNetRelease_t = ...
        -leafStorageRate_t;


    leafRelease_t = ...
        max( ...
        leafNetRelease_t, ...
        0);


    leafRecharge_t = ...
        max( ...
        -leafNetRelease_t, ...
        0);


    JxylemToLeaf(t) = ...
        JxylemToLeaf_t;


    leafHydraulicStorageAmount(t) = ...
        leafStorageAmount_t;


    leafHydraulicStorageRate(t) = ...
        leafStorageRate_t;


    leafStorageNetReleaseFlux(t) = ...
        leafNetRelease_t;


    leafStorageReleaseFlux(t) = ...
        leafRelease_t;


    leafStorageRechargeFlux(t) = ...
        leafRecharge_t;


    leafCapacitanceEffective(t) = ...
        leafCapacitanceEffective_t;


    leafStorageBoundActive(t) = ...
        leafStorageBound_t;

    %% ------------------------------------------------------------
    % Accepted leaf pressure-volume diagnostics
    %
    % Re-evaluate the constitutive relation using the accepted
    % psiLeaf solution. This ensures diagnostics are populated at
    % every timestep, including t = 1.
    % ------------------------------------------------------------

    if ...
            leafCapacitanceEnabled ...
            && any(leafStorageModel == [ ...
            "pressurevolume", ...
            "pressure-volume", ...
            "pv"])

        SleafAccepted = ...
            leaf_storage_state( ...
            psiL, ...
            leafStorageCfg);


        leafRWC(t) = ...
            SleafAccepted.relativeWaterContent;


        leafSymplasticRWC(t) = ...
            SleafAccepted.symplasticRelativeWaterContent;


        leafOsmoticPsi(t) = ...
            SleafAccepted.osmoticPotential_MPa;


        leafTurgorPsi(t) = ...
            SleafAccepted.turgorPotential_MPa;


        leafApoplasticWater(t) = ...
            SleafAccepted.apoplasticWaterAmount_mol_m2;


        leafSymplasticWater(t) = ...
            SleafAccepted.symplasticWaterAmount_mol_m2;

    end

    if ...
            t == 1 ...
            && leafCapacitanceEnabled

        leafHydraulicStorageInitialActual = ...
            leafStorageAmount_t;

    end


    %% ------------------------------------------------------------
    % Accepted hydraulic-storage state
    % ------------------------------------------------------------

    if ...
            storageModel == "twocomponentpiecewise" ...
            || storageModel == "twocomponenthysteretic"

        hydraulicStorageAmount(t) = ...
            state.storageAmount;


        hydraulicStorageLumenAmount(t) = ...
            state.storageLumenAmount;


        hydraulicStorageLivingAmount(t) = ...
            state.storageLivingAmount;


        xylemCapacitanceEffective(t) = ...
            state.capacitanceEffective;


        xylemCapacitanceLumen(t) = ...
            state.capacitanceLumen;


        xylemCapacitanceLiving(t) = ...
            state.capacitanceLiving;

    end

    %% ------------------------------------------------------------
    % Layer-specific uptake
    % ------------------------------------------------------------

    if isempty(layer)
        %do nothing
    else
        Ji = layer.flux;
    end

    if any(~isfinite(Ji))

        solverFallback(t) = true;

        break

    end


    totalJ = ...
        sum(Ji);


    Juptake = ...
        max(Ji,0);


    Jrelease = ...
        max(-Ji,0);


    totalUptake = ...
        sum(Juptake);


    totalRelease = ...
        sum(Jrelease);


    if totalUptake > 0

        wi = ...
            Juptake ...
            ./ totalUptake;

    else

        wi = ...
            NaN(nLayer,1);

    end


    %% ------------------------------------------------------------
    % Initialize storage reference state
    % ------------------------------------------------------------
    

    if t == 1

        psiX0 = ...
            psiX;


        if storageModel == "constantcapacitance"

            hydraulicStorageInitialActual = ...
                hydraulicStorageInitial;


            hydraulicStorageAmount(1) = ...
                hydraulicStorageInitialActual;


        elseif ...
                storageModel == "twocomponentpiecewise" ...
                || storageModel == "twocomponenthysteretic"

            % Anatomical storage fields define maximum capacity.
            % Actual initial storage depends on initial PsiX.

            hydraulicStorageInitialActual = ...
                state.storageAmount;

        end

    end


    %% ------------------------------------------------------------
    % Capacitance fluxes
    %
    % storageRate:
    %
    %   dS/dt = Jroot - E
    %
    % netStorageRelease:
    %
    %   E - Jroot
    %
    % Positive net release = storage supplying transpiration.
    % ------------------------------------------------------------

    storageRate_t = ...
        state.storageRate;

    %% ------------------------------------------------------------
    % Component-specific storage rates
    %
    % Positive rate  = storage recharge
    % Negative rate  = storage release
    % ------------------------------------------------------------

    if ...
            storageModel == "twocomponentpiecewise" ...
            || storageModel == "twocomponenthysteretic"

        if t == 1

            storageLumenRate_t = ...
                0;

            storageLivingRate_t = ...
                0;

        else

            storageLumenRate_t = ...
                ( ...
                hydraulicStorageLumenAmount(t) ...
                - hydraulicStorageLumenAmount(t-1) ...
                ) ...
                ./ dtSeconds;


            storageLivingRate_t = ...
                ( ...
                hydraulicStorageLivingAmount(t) ...
                - hydraulicStorageLivingAmount(t-1) ...
                ) ...
                ./ dtSeconds;

        end


        %% Positive release fluxes

        lumenRelease_t = ...
            max( ...
            -storageLumenRate_t, ...
            0);


        livingRelease_t = ...
            max( ...
            -storageLivingRate_t, ...
            0);


        %% Positive recharge fluxes

        lumenRecharge_t = ...
            max( ...
            storageLumenRate_t, ...
            0);


        livingRecharge_t = ...
            max( ...
            storageLivingRate_t, ...
            0);


    else

        storageLumenRate_t = ...
            NaN;

        storageLivingRate_t = ...
            NaN;


        lumenRelease_t = ...
            NaN;

        lumenRecharge_t = ...
            NaN;


        livingRelease_t = ...
            NaN;

        livingRecharge_t = ...
            NaN;

    end


    netRelease_t = ...
        -storageRate_t;


    release_t = ...
        max(netRelease_t,0);


    recharge_t = ...
        max(-netRelease_t,0);

    %% ------------------------------------------------------------
    % Component-storage closure
    % ------------------------------------------------------------

    if ...
            storageModel == "twocomponentpiecewise" ...
            || storageModel == "twocomponenthysteretic"

        componentStorageRateResidual_t = ...
            storageRate_t ...
            - ( ...
            storageLumenRate_t ...
            + storageLivingRate_t);


        componentStorageNetFluxResidual_t = ...
            netRelease_t ...
            - ( ...
            (lumenRelease_t - lumenRecharge_t) ...
            + ...
            (livingRelease_t - livingRecharge_t));

    else

        componentStorageRateResidual_t = ...
            NaN;

        componentStorageNetFluxResidual_t = ...
            NaN;

    end


    if E_t > 0

        storageFraction_t = ...
            release_t ./ E_t;

    else

        storageFraction_t = NaN;

    end


    %% ------------------------------------------------------------
    % Store primary states
    % ------------------------------------------------------------

    psiRoot(t) = ...
        psiR;

    psiLeaf(t) = ...
        psiL;

    psiXylem(t) = ...
        psiX;

    E(t) = ...
        E_t;

    gs(t) = ...
        gs_t;

    rootFlux(:,t) = ...
        Ji;

    rootWeights(:,t) = ...
        wi;

    KstemOut(t) = ...
        Kstem_t;

    KleafOut(t) = ...
        Kleaf_t;


    %% ------------------------------------------------------------
    % Store capacitance states
    % ------------------------------------------------------------

    requiredStemFlow(t) = ...
        state.requiredStemFlow;

    actualStemFlow(t) = ...
        state.actualStemFlow;

    storageRate(t) = ...
        storageRate_t;

    netStorageReleaseFlux(t) = ...
        netRelease_t;

    storageReleaseFlux(t) = ...
        release_t;

    storageRechargeFlux(t) = ...
        recharge_t;

    storageContributionFraction(t) = ...
        storageFraction_t;

    storageLumenRate(t) = ...
        storageLumenRate_t;


    storageLivingRate(t) = ...
        storageLivingRate_t;


    storageLumenReleaseFlux(t) = ...
        lumenRelease_t;


    storageLumenRechargeFlux(t) = ...
        lumenRecharge_t;


    storageLivingReleaseFlux(t) = ...
        livingRelease_t;


    storageLivingRechargeFlux(t) = ...
        livingRecharge_t;


    componentStorageRateResidual(t) = ...
        componentStorageRateResidual_t;


    componentStorageNetFluxResidual(t) = ...
        componentStorageNetFluxResidual_t;


    if capacitanceEnabled

        switch storageModel


            case "constantcapacitance"

                deltaXylemStorage(t) = ...
                    Cx ...
                    .* (psiX - psiX0);


                hydraulicStorageAmount(t) = ...
                    hydraulicStorageInitialActual ...
                    + deltaXylemStorage(t);


            case { ...
                    "twocomponentpiecewise", ...
                    "twocomponenthysteretic"}

                deltaXylemStorage(t) = ...
                    hydraulicStorageAmount(t) ...
                    - hydraulicStorageInitialActual;

        end

        if ...
                storageModel == "twocomponentpiecewise" ...
                || storageModel == "twocomponenthysteretic"

            storageLumenBranchCode(t) = ...
                state.storageLumenBranchCode;

        end


        xylemStorageWaterDeficit(t) = ...
            -deltaXylemStorage(t);


    else

        deltaXylemStorage(t) = ...
            0;

        xylemStorageWaterDeficit(t) = ...
            0;

    end


    %% ------------------------------------------------------------
    % Nonlinear soil-root diagnostics
    % ------------------------------------------------------------

    psiRootSurfaceOut(:,t) = ...
        layer.psiRootSurface;

    KrootTissueOut(:,t) = ...
        layer.KrootTissue;

    KsoilPathOut(:,t) = ...
        layer.KsoilEquivalent;

    KeffectiveOut(:,t) = ...
        layer.Keffective;

    soilLimitOut(:,t) = ...
        layer.soilLimitationFraction;

    rootLimitOut(:,t) = ...
        layer.rootLimitationFraction;

    soilKrelOut(:,t) = ...
        layer.relativeSoilConductivity;

    KsoilUnsatOut(:,t) = ...
        layer.KsoilUnsat_m_s;

    KsoilPathBulkOut(:,t) = ...
        layer.KsoilPathBulk_m_s;

    soilFluxResidualOut(:,t) = ...
        layer.soilFluxResidual;


    %% ------------------------------------------------------------
    % Hydraulic mass balance
    %
    % total root uptake + storage release = transpiration
    % ------------------------------------------------------------

    massBalanceResidual(t) = ...
        totalJ ...
        + netRelease_t ...
        + leafNetRelease_t ...
        - E_t;


    %% Stem-flow closure

    stemFlowResidual(t) = ...
        state.actualStemFlow ...
        - state.requiredStemFlow;


    % Xylem-node closure:
    %
    % Jroot - JxylemToLeaf - dSx/dt = 0

    capacitanceMassBalanceResidual(t) = ...
        totalJ ...
        - JxylemToLeaf_t ...
        - storageRate_t;

    rootUptakeFlux(:,t) = ...
        Juptake;


    rootReleaseFlux(:,t) = ...
        Jrelease;


    totalRootUptakeFlux(t) = ...
        totalUptake;


    totalRootReleaseFlux(t) = ...
        totalRelease;


    hydraulicallyFeasible(t) = true;

end


%% ========================================================================
% Package outputs
% ========================================================================

hydro.leafRelativeWaterContent = ...
    leafRWC;


hydro.leafSymplasticRelativeWaterContent = ...
    leafSymplasticRWC;


hydro.leafOsmoticPotential_MPa = ...
    leafOsmoticPsi;


hydro.leafTurgorPotential_MPa = ...
    leafTurgorPsi;


hydro.leafApoplasticWaterAmount_mol_m2 = ...
    leafApoplasticWater;


hydro.leafSymplasticWaterAmount_mol_m2 = ...
    leafSymplasticWater;

hydro.psiRoot = ...
    psiRoot;

hydro.psiLeaf = ...
    psiLeaf;

%% ========================================================================
% Two-domain leaf hydraulic states
% ========================================================================

hydro.leafParallelHydraulicsEnabled = ...
    logical(leafNetCfg.enabled);


hydro.psiSymplast = ...
    psiSymplast;


hydro.psiApoplast = ...
    psiApoplast;


hydro.psiEvaporationSite = ...
    psiEvaporationSite;


hydro.JxylemToApoplast = ...
    JxylemToApoplast;


hydro.JapoplastToSymplast = ...
    JapoplastToSymplast;


hydro.JsymplastToEvaporation = ...
    JsymplastToEvaporation;


hydro.JapoplastToEvaporation = ...
    JapoplastToEvaporation;


hydro.symplasticAdvectiveFraction = ...
    symplasticAdvectiveFractionHydraulic;


hydro.apoplasticAdvectiveFraction = ...
    apoplasticAdvectiveFractionHydraulic;


hydro.KxylemApoplast = ...
    KxylemApoplastOut;


hydro.leafParallelApoplastClosure_mol_m2_s = ...
    leafParallelApoplastClosure;


hydro.leafParallelSymplastClosure_mol_m2_s = ...
    leafParallelSymplastClosure;


hydro.leafParallelEvaporationClosure_mol_m2_s = ...
    leafParallelEvaporationClosure;


hydro.leafParallelWholeClosure_mol_m2_s = ...
    leafParallelWholeClosure;


hydro.leafParallelStorageRateResidual_mol_m2_s = ...
    leafParallelStorageRateResidual;

hydro.psiXylem = ...
    psiXylem;


hydro.E = ...
    E;

hydro.gs = ...
    gs;


hydro.rootFlux = ...
    rootFlux;

hydro.rootWeights = ...
    rootWeights;


%% ------------------------------------------------------------------------
% Root uptake
% -------------------------------------------------------------------------

hydro.rootUptakeTotal = ...
    NaN(1,nTime);


validFlux = ...
    all(isfinite(rootFlux),1);


hydro.rootUptakeTotal(validFlux) = ...
    sum( ...
        rootFlux(:,validFlux), ...
        1);

%% ------------------------------------------------------------------------
% Explicit water-flux bookkeeping for xylem isotope coupling
% -------------------------------------------------------------------------

hydro.Jroot = ...
    hydro.rootUptakeTotal;

% Net atmospheric water loss

hydro.Jtranspiration = ...
    E;


% Water delivered from terminal xylem into the leaf

hydro.JxylemToLeaf = ...
    JxylemToLeaf;



% Jleaf now represents actual xylem -> leaf delivery.

hydro.Jleaf = ...
    JxylemToLeaf;

% Positive = storage -> mobile xylem
hydro.JstorageRelease = ...
    storageReleaseFlux;

% Positive = mobile xylem -> storage
hydro.JstorageRecharge = ...
    storageRechargeFlux;

% Signed net storage contribution:
%
% positive = release
% negative = recharge

hydro.JstorageNet = ...
    storageReleaseFlux ...
    - storageRechargeFlux;

% Hydraulic flow identity:
%
% Jroot + JstorageNet = Jleaf

hydro.storageFluxClosureResidual = ...
    hydro.Jroot ...
    + hydro.JstorageNet ...
    - hydro.Jleaf;

%% ------------------------------------------------------------------------
% Leaf hydraulic storage fluxes
% -------------------------------------------------------------------------

hydro.JleafStorageRelease = ...
    leafStorageReleaseFlux;


hydro.JleafStorageRecharge = ...
    leafStorageRechargeFlux;


hydro.JleafStorageNet = ...
    leafStorageReleaseFlux ...
    - leafStorageRechargeFlux;


hydro.leafStorageFluxClosureResidual = ...
    hydro.JxylemToLeaf ...
    + hydro.JleafStorageNet ...
    - hydro.E;


%% Whole-plant identity:
%
% Jroot
% + xylem storage release
% + leaf storage release
% = transpiration

hydro.wholePlantStorageFluxClosureResidual = ...
    hydro.Jroot ...
    + hydro.JstorageNet ...
    + hydro.JleafStorageNet ...
    - hydro.E;

hydro.leafCapacitanceEnabled = ...
    leafCapacitanceEnabled;

hydro.diagnosticLeafStorageFluxDisabled = ...
    diagnosticDisableLeafStorageFlux;


hydro.leafHydraulicStorageInitial_mol_m2 = ...
    leafHydraulicStorageInitialActual;


hydro.leafHydraulicStorageAmount_mol_m2 = ...
    leafHydraulicStorageAmount;


hydro.leafHydraulicStorageRate_mol_m2_s = ...
    leafHydraulicStorageRate;


hydro.leafCapacitanceEffective_mol_m2_MPa = ...
    leafCapacitanceEffective;


hydro.leafStorageBoundActive = ...
    leafStorageBoundActive;


hydro.leafCapacitanceMassBalanceResidual = ...
    leafStorageFluxClosureResidual;

%% ------------------------------------------------------------------------
% Absolute hydraulic-storage state
% -------------------------------------------------------------------------

if capacitanceEnabled

    hydro.hydraulicStorageInitial_mol_m2 = ...
        hydraulicStorageInitialActual;


    hydro.hydraulicStorageAmount_mol_m2 = ...
        hydraulicStorageAmount;


    hydro.hydraulicStorageFractionRemaining = ...
    hydraulicStorageAmount ...
    ./ hydraulicStorageInitialActual;


    hydro.hydraulicStorageCapacityViolated = ...
        hydraulicStorageAmount <= 0;


    hydro.hydraulicStorageCapacity_mol_m2 = ...
        hydraulicStorageCapacity;


    if isfinite(hydraulicStorageCapacity)

        hydro.hydraulicStorageFractionOfCapacity = ...
            hydraulicStorageAmount ...
            ./ hydraulicStorageCapacity;

    else

        hydro.hydraulicStorageFractionOfCapacity = ...
            NaN(1,nTime);

    end


    if ...
            storageModel == "twocomponentpiecewise" ...
            || storageModel == "twocomponenthysteretic"

        hydro.hydraulicStorageLumenAmount_mol_m2 = ...
            hydraulicStorageLumenAmount;


        hydro.hydraulicStorageLivingAmount_mol_m2 = ...
            hydraulicStorageLivingAmount;

    else

        hydro.hydraulicStorageLumenAmount_mol_m2 = ...
            NaN(1,nTime);


        hydro.hydraulicStorageLivingAmount_mol_m2 = ...
            NaN(1,nTime);

    end

else

    hydro.hydraulicStorageInitial_mol_m2 = ...
        NaN;


    hydro.hydraulicStorageAmount_mol_m2 = ...
        NaN(1,nTime);


    hydro.hydraulicStorageFractionRemaining = ...
        NaN(1,nTime);


    hydro.hydraulicStorageCapacityViolated = ...
        false(1,nTime);


    hydro.hydraulicStorageCapacity_mol_m2 = ...
        NaN;

end

%% ------------------------------------------------------------------------
% Capacitance outputs
% -------------------------------------------------------------------------

hydro.xylemCapacitanceEnabled = ...
    capacitanceEnabled;

hydro.diagnosticXylemStorageFluxDisabled = ...
    diagnosticDisableXylemStorageFlux;


hydro.xylemCapacitance_mol_m2_MPa = ...
    Cx;

hydro.xylemCapacitanceIsStateDependent = ...
    capacitanceEnabled ...
    && any(storageModel == [ ...
        "twocomponentpiecewise", ...
        "twocomponenthysteretic"]);

hydro.storageModel = ...
    char(storageModel);


if ...
        storageModel == "twocomponentpiecewise" ...
        || storageModel == "twocomponenthysteretic"

    hydro.xylemCapacitanceEffective_mol_m2_MPa = ...
        xylemCapacitanceEffective;


    hydro.xylemCapacitanceLumen_mol_m2_MPa = ...
        xylemCapacitanceLumen;


    hydro.xylemCapacitanceLiving_mol_m2_MPa = ...
        xylemCapacitanceLiving;


else

    hydro.xylemCapacitanceEffective_mol_m2_MPa = ...
        repmat(Cx,1,nTime);


    hydro.xylemCapacitanceLumen_mol_m2_MPa = ...
        NaN(1,nTime);


    hydro.xylemCapacitanceLiving_mol_m2_MPa = ...
        NaN(1,nTime);

end

hydro.capacitanceParameterization = ...
    capacitanceParameterization;


hydro.requiredStemFlow = ...
    requiredStemFlow;


hydro.actualStemFlow = ...
    actualStemFlow;


% Positive = recharge / increasing stored water

hydro.hydraulicStorageRate_mol_m2_s = ...
    storageRate;


% Positive = storage supplying transpiration

hydro.hydraulicStorageNetReleaseFlux_mol_m2_s = ...
    netStorageReleaseFlux;


hydro.hydraulicStorageReleaseFlux_mol_m2_s = ...
    storageReleaseFlux;


hydro.hydraulicStorageRechargeFlux_mol_m2_s = ...
    storageRechargeFlux;


hydro.hydraulicStorageContributionFraction = ...
    storageContributionFraction;


hydro.deltaHydraulicStorage_mol_m2 = ...
    deltaXylemStorage;


hydro.hydraulicStorageWaterDeficit_mol_m2 = ...
    xylemStorageWaterDeficit;


hydro.capacitanceMassBalanceResidual = ...
    capacitanceMassBalanceResidual;


hydro.stemFlowResidual = ...
    stemFlowResidual;


%% ------------------------------------------------------------------------
% Dynamic capacitance timescale
% -------------------------------------------------------------------------

hydro.xylemHydraulicTimeConstantDays = ...
    NaN(1,nTime);


validStemK = ...
    isfinite(KstemOut) ...
    & KstemOut > 0;


if capacitanceEnabled

    if storageModel == "constantcapacitance"

        hydro.xylemHydraulicTimeConstantDays(validStemK) = ...
            Cx ...
            ./ KstemOut(validStemK) ...
            ./ secondsPerDay;


    elseif ...
            storageModel == "twocomponentpiecewise" ...
            || storageModel == "twocomponenthysteretic"

        validTau = ...
            validStemK ...
            & isfinite(xylemCapacitanceEffective);


        hydro.xylemHydraulicTimeConstantDays(validTau) = ...
            xylemCapacitanceEffective(validTau) ...
            ./ KstemOut(validTau) ...
            ./ secondsPerDay;

    end

end


%% ------------------------------------------------------------------------
% Cumulative storage release/recharge
% -------------------------------------------------------------------------

hydro.cumulativeStorageRelease_mol_m2 = ...
    NaN(1,nTime);

hydro.cumulativeStorageRecharge_mol_m2 = ...
    NaN(1,nTime);


validStorage = ...
    isfinite(storageReleaseFlux) ...
    & isfinite(storageRechargeFlux);


if any(validStorage)

    lastValidStorage = ...
        find(validStorage,1,'last');


    ii = ...
        1:lastValidStorage;


    %% ------------------------------------------------------------
    % Singleton-safe cumulative storage bookkeeping
    %
    % With only one valid timestep there is no elapsed interval to
    % integrate, so cumulative release/recharge are exactly zero.
    % ------------------------------------------------------------

    if numel(ii) == 1

        hydro.cumulativeStorageRelease_mol_m2(ii) = ...
            0;


        hydro.cumulativeStorageRecharge_mol_m2(ii) = ...
            0;


    else

        timeSeconds = ...
            timeDays(ii) ...
            .* secondsPerDay;


        hydro.cumulativeStorageRelease_mol_m2(ii) = ...
            cumtrapz( ...
                timeSeconds, ...
                storageReleaseFlux(ii), ...
                2);


        hydro.cumulativeStorageRecharge_mol_m2(ii) = ...
            cumtrapz( ...
                timeSeconds, ...
                storageRechargeFlux(ii), ...
                2);

    end

end


%% ------------------------------------------------------------------------
% Whole-plant mass balance
% -------------------------------------------------------------------------

hydro.massBalanceResidual = ...
    massBalanceResidual;


%% ------------------------------------------------------------------------
% Conductance and vulnerability
% -------------------------------------------------------------------------

hydro.Kstem = ...
    KstemOut;

hydro.Kleaf = ...
    KleafOut;


hydro.PLCstem = ...
    100 ...
    .* (1 - KstemOut ./ plantState.KstemMax);


hydro.PLCleaf = ...
    100 ...
    .* (1 - KleafOut ./ plantState.KleafMax);


%% ------------------------------------------------------------------------
% Solver diagnostics
% -------------------------------------------------------------------------

hydro.hydraulicResidual = ...
    hydraulicResidualOut;

hydro.solverFallback = ...
    solverFallback;

hydro.nHydraulicRoots = ...
    nHydraulicRoots;

hydro.hydraulicallyFeasible = ...
    hydraulicallyFeasible;

hydro.allHydraulicRoots = ...
    allHydraulicRoots;

hydro.failureDiagnostics = ...
    failureDiagnostics;

hydro.branchGs = ...
    branchGs;

hydro.branchE = ...
    branchE;

hydro.branchKleaf = ...
    branchKleaf;


%% ------------------------------------------------------------------------
% Nonlinear soil-root pathway
% -------------------------------------------------------------------------

hydro.psiRootSurface = ...
    psiRootSurfaceOut;


hydro.KrootTissue = ...
    KrootTissueOut;


hydro.KsoilPath = ...
    KsoilPathOut;


hydro.KrootEffective = ...
    KeffectiveOut;


hydro.soilLimitationFraction = ...
    soilLimitOut;


hydro.rootLimitationFraction = ...
    rootLimitOut;


hydro.relativeSoilConductivity = ...
    soilKrelOut;


hydro.KsoilUnsat_m_s = ...
    KsoilUnsatOut;


hydro.KsoilPathBulk_m_s = ...
    KsoilPathBulkOut;


hydro.soilFluxResidual = ...
    soilFluxResidualOut;


%% ------------------------------------------------------------------------
% Static root geometry
% -------------------------------------------------------------------------

hydro.rootRadius_m = ...
    repmat( ...
        rootPath.rootRadius_m, ...
        1,nTime);


hydro.rOuter_m = ...
    repmat( ...
        rootPath.rOuter_m, ...
        1,nTime);


hydro.logGeometry = ...
    repmat( ...
        rootPath.logGeometry, ...
        1,nTime);


hydro.rootLengthPerLeafArea = ...
    repmat( ...
        rootPath.rootLengthPerLeafArea, ...
        1,nTime);


hydro.rootSurfaceAreaPerLeafArea = ...
    repmat( ...
        rootPath.rootSurfaceAreaPerLeafArea, ...
        1,nTime);


hydro.geometryRootFraction = ...
    repmat( ...
        rootPath.geometryRootFraction, ...
        1,nTime);


%% ------------------------------------------------------------------------
% Potential drops
% -------------------------------------------------------------------------

hydro.soilPotentialDrop = ...
    psiSoil ...
    - psiRootSurfaceOut;


hydro.rootTissuePotentialDrop = ...
    psiRootSurfaceOut ...
    - repmat(psiRoot,nLayer,1);


hydro.stemPotentialDrop = ...
    psiRoot ...
    - psiXylem ...
    - gravityMPa;


hydro.leafPotentialDrop = ...
    psiXylem ...
    - psiLeaf;

hydro.xylemApoplastPotentialDrop = ...
    psiXylem ...
    - psiApoplast;


hydro.apoplastSymplastPotentialDrop = ...
    psiApoplast ...
    - psiSymplast;


hydro.symplastEvaporationPotentialDrop = ...
    psiSymplast ...
    - psiEvaporationSite;


hydro.apoplastEvaporationPotentialDrop = ...
    psiApoplast ...
    - psiEvaporationSite;


hydro.rootPathwayMethod = ...
    'matric-flux-potential';


%% ------------------------------------------------------------------------
% Effective root-system conductance
% -------------------------------------------------------------------------

hydro.KrootSystem = ...
    NaN(1,nTime);


validRootK = ...
    any(isfinite(KeffectiveOut),1);


hydro.KrootSystem(validRootK) = ...
    sum( ...
        KeffectiveOut(:,validRootK), ...
        1, ...
        'omitnan');


hydro.RrootSystem = ...
    NaN(1,nTime);


validRootSystem = ...
    isfinite(hydro.KrootSystem) ...
    & hydro.KrootSystem > 0;


hydro.RrootSystem(validRootSystem) = ...
    1 ./ hydro.KrootSystem(validRootSystem);


%% ------------------------------------------------------------------------
% Whole-plant resistance diagnostics
%
% These remain instantaneous resistance diagnostics.
% Capacitance is not itself a resistance.
% -------------------------------------------------------------------------

hydro.Rstem = ...
    NaN(1,nTime);

hydro.Rleaf = ...
    NaN(1,nTime);


validStem = ...
    isfinite(hydro.Kstem) ...
    & hydro.Kstem > 0;


validLeaf = ...
    isfinite(hydro.Kleaf) ...
    & hydro.Kleaf > 0;


hydro.Rstem(validStem) = ...
    1 ./ hydro.Kstem(validStem);


hydro.Rleaf(validLeaf) = ...
    1 ./ hydro.Kleaf(validLeaf);


hydro.Rtotal = ...
    NaN(1,nTime);


validR = ...
    isfinite(hydro.RrootSystem) ...
    & isfinite(hydro.Rstem) ...
    & isfinite(hydro.Rleaf);


hydro.Rtotal(validR) = ...
    hydro.RrootSystem(validR) ...
    + hydro.Rstem(validR) ...
    + hydro.Rleaf(validR);


hydro.rootResistanceFraction = ...
    NaN(1,nTime);

hydro.stemResistanceFraction = ...
    NaN(1,nTime);

hydro.leafResistanceFraction = ...
    NaN(1,nTime);


hydro.rootResistanceFraction(validR) = ...
    hydro.RrootSystem(validR) ...
    ./ hydro.Rtotal(validR);


hydro.stemResistanceFraction(validR) = ...
    hydro.Rstem(validR) ...
    ./ hydro.Rtotal(validR);


hydro.leafResistanceFraction(validR) = ...
    hydro.Rleaf(validR) ...
    ./ hydro.Rtotal(validR);


%% ------------------------------------------------------------------------
% Root-source accessibility
% -------------------------------------------------------------------------

hydro.activeRootLayers = ...
    isfinite(rootFlux) ...
    & rootFlux > 0;


hydro.nActiveRootLayers = ...
    NaN(1,nTime);


validFluxTime = ...
    any(isfinite(rootFlux),1);


hydro.nActiveRootLayers(validFluxTime) = ...
    sum( ...
        hydro.activeRootLayers(:,validFluxTime), ...
        1);


hydro.rootAccessMargin = ...
    soil.psiMPa ...
    - repmat(hydro.psiRoot,nLayer,1);


hydro.rootPath = ...
    rootPath;

hydro.storageLumenBranchCode = ...
    storageLumenBranchCode;


%% ------------------------------------------------------------------------
% Component-specific hydraulic storage fluxes
% -------------------------------------------------------------------------

hydro.hydraulicStorageLumenRate_mol_m2_s = ...
    storageLumenRate;


hydro.hydraulicStorageLivingRate_mol_m2_s = ...
    storageLivingRate;


% Positive component release:
% storage -> mobile xylem

hydro.JstorageLumenRelease = ...
    storageLumenReleaseFlux;


hydro.JstorageLivingRelease = ...
    storageLivingReleaseFlux;


% Positive component recharge:
% mobile xylem -> storage

hydro.JstorageLumenRecharge = ...
    storageLumenRechargeFlux;


hydro.JstorageLivingRecharge = ...
    storageLivingRechargeFlux;


% Signed component contribution to mobile xylem:
% positive = release
% negative = recharge

hydro.JstorageLumenNet = ...
    storageLumenReleaseFlux ...
    - storageLumenRechargeFlux;


hydro.JstorageLivingNet = ...
    storageLivingReleaseFlux ...
    - storageLivingRechargeFlux;


hydro.componentStorageRateResidual = ...
    componentStorageRateResidual;


hydro.componentStorageNetFluxResidual = ...
    componentStorageNetFluxResidual;

hydro.rootFlux = ...
    rootFlux;


hydro.rootFluxGround_mol_m2_s = ...
    plantState.leafAreaIndex ...
    .* rootFlux;


hydro.rootFluxAreaBasis = ...
    'leaf area';


hydro.rootFluxGroundAreaBasis = ...
    'ground area';


hydro.leafAreaIndex = ...
    plantState.leafAreaIndex;


hydro.transpirationAreaBasis = ...
    'leaf area';

hydro.usedPreviousStateContinuation = ...
    usedContinuation;


hydro.previousStateResidual = ...
    continuationResidual;


hydro.previousStateValid = ...
    continuationStateValid;

end


%% =========================================================================
function [residual,state] = hydraulicResidualCapacitance( ...
    psiLeaf, ...
    psiSoil, ...
    rootPath, ...
    met, ...
    plantState, ...
    t, ...
    psiXPrevious, ...
    storagePrevious, ...
    storageLumenPrevious, ...
    dtSeconds, ...
    Cx, ...
    useCapacitance, ...
    storageModel, ...
    xylemTraits, ...
    leafStoragePrevious, ...
    useLeafCapacitance, ...
    leafStorageCfg, ...
    leafNetCfg, ...
    leafSymplasticWaterPrevious, ...
    gravityMPa, ...
    leafTraits, ...
    rootDemandBracketCache, ...
    psiRootBranchReference)
%HYDRAULICRESIDUALCAPACITANCE
%
% Hydraulic topology:
%
% soil -> root -> stem -> xylem storage -> leaf storage -> atmosphere
%
%
% Leaf node:
%
%   JxylemToLeaf =
%       E + dSleaf/dt
%
%
% Leaf conductance:
%
%   JxylemToLeaf =
%       Kleaf * (psiX - psiLeaf)
%
%
% Xylem node:
%
%   Jroot =
%       JxylemToLeaf + dSxylem/dt
%
%
% Outer residual:
%
%   Jstem - JrootRequired


%% ========================================================================
% Default invalid state
% ========================================================================

state.valid = ...
    false;

% Diagnostic classification for invalid candidate states. This field does
% not alter hydraulic feasibility; it only records the first rejected
% constraint encountered by the existing equations.
state.rejectCode = ...
    "unclassified_invalid";


state.E = ...
    NaN;


state.gs = ...
    NaN;


state.psiRoot = ...
    NaN;


state.psiXylem = ...
    NaN;


state.Kstem = ...
    NaN;


state.Kleaf = ...
    NaN;


%% Xylem storage

state.storageRate = ...
    NaN;


state.storageAmount = ...
    NaN;


state.storageLumenAmount = ...
    NaN;


state.storageLivingAmount = ...
    NaN;


state.capacitanceEffective = ...
    NaN;


state.capacitanceLumen = ...
    NaN;


state.capacitanceLiving = ...
    NaN;


state.storageLumenBranchCode = ...
    NaN;


state.storageLumenReleaseActive = ...
    false;


state.storageLumenRefillActive = ...
    false;


%% Leaf storage

state.leafStorageRate = ...
    NaN;


state.leafStorageAmount = ...
    NaN;


state.leafCapacitanceEffective = ...
    NaN;


state.leafStorageBoundActive = ...
    false;

%% Leaf pressure-volume diagnostics

state.leafRelativeWaterContent = ...
    NaN;


state.leafSymplasticRelativeWaterContent = ...
    NaN;


state.leafOsmoticPotential_MPa = ...
    NaN;


state.leafTurgorPotential_MPa = ...
    NaN;


state.leafApoplasticWaterAmount_mol_m2 = ...
    NaN;


state.leafSymplasticWaterAmount_mol_m2 = ...
    NaN;

state.JxylemToLeaf = ...
    NaN;


%% Flow diagnostics

state.requiredStemFlow = ...
    NaN;


state.actualStemFlow = ...
    NaN;


state.layer = [];


residual = ...
    NaN;

%% Parallel leaf network

state.psiSymplast = ...
    NaN;


state.psiApoplast = ...
    NaN;


state.psiEvaporationSite = ...
    NaN;


state.JxylemToApoplast = ...
    NaN;


state.JapoplastToSymplast = ...
    NaN;


state.JsymplastToEvaporation = ...
    NaN;


state.JapoplastToEvaporation = ...
    NaN;


state.symplasticAdvectiveFraction = ...
    NaN;


state.apoplasticAdvectiveFraction = ...
    NaN;


state.leafParallelApoplastClosure = ...
    NaN;


state.leafParallelSymplastClosure = ...
    NaN;


state.leafParallelEvaporationClosure = ...
    NaN;


state.leafParallelWholeClosure = ...
    NaN;


state.leafParallelStorageRateResidual = ...
    NaN;


state.KxylemApoplast = ...
    NaN;


%% ========================================================================
% Candidate leaf potential
% ========================================================================

if ~isfinite(psiLeaf)

    state.rejectCode = "nonfinite_psi_leaf";
    return

end


%% ========================================================================
% Atmospheric demand
% ========================================================================

[E,gs] = ...
    transpirationDemand( ...
        psiLeaf, ...
        met, ...
        plantState, ...
        t);


if ...
        ~isfinite(E) ...
        || E < 0 ...
        || ~isfinite(gs)

    state.rejectCode = "invalid_atmospheric_demand";
    return

end


%% ========================================================================
% Dynamic leaf conductance
% ========================================================================

%% ========================================================================
% Xylem -> apoplast conductance
%
% In legacy mode this reproduces the previously validated
% vulnerability-dependent Kleaf exactly.
% ========================================================================

if leafNetCfg.enabled

    KXA = ...
        resolveLeafXylemApoplastConductance( ...
            psiLeaf, ...
            plantState, ...
            leafNetCfg);


    % Compatibility:
    % state.Kleaf now represents KXA when parallel hydraulics are active.

    Kleaf = ...
        KXA;

else

    Kleaf = ...
        hydraulic_vulnerability( ...
            psiLeaf, ...
            plantState.KleafMax, ...
            plantState.P50Leaf, ...
            plantState.vulnerabilityShapeLeaf);


    KXA = ...
        NaN;

end


if ...
        ~isfinite(Kleaf) ...
        || Kleaf <= 0

    state.rejectCode = "invalid_leaf_conductance";
    return

end


if ...
        ~isfinite(Kleaf) ...
        || Kleaf <= 0

    state.rejectCode = "invalid_leaf_conductance";
    return

end


%% ========================================================================
% Leaf hydraulic storage
% ========================================================================

if isfield(leafStorageCfg,'enabled') ...
        && logical(leafStorageCfg.enabled)

    Sleaf = ...
        leaf_storage_state( ...
        psiLeaf, ...
        leafStorageCfg);


    leafStorageAmountCandidate = ...
        Sleaf.amount_mol_m2;


    leafCapacitanceCandidate = ...
        Sleaf.capacitanceEffective_mol_m2_MPa;


    leafStorageBoundCandidate = ...
        Sleaf.boundActive;

    %% ------------------------------------------------------------
    % Optional pressure-volume diagnostics
    % ------------------------------------------------------------

    leafRWCCandidate = ...
        NaN;


    leafSymplasticRWCCandidate = ...
        NaN;


    leafOsmoticPsiCandidate = ...
        NaN;


    leafTurgorPsiCandidate = ...
        NaN;


    leafApoplasticWaterCandidate = ...
        NaN;


    leafSymplasticWaterCandidate = ...
        NaN;


    if isfield(Sleaf,'relativeWaterContent')

        leafRWCCandidate = ...
            Sleaf.relativeWaterContent;


        leafSymplasticRWCCandidate = ...
            Sleaf.symplasticRelativeWaterContent;


        leafOsmoticPsiCandidate = ...
            Sleaf.osmoticPotential_MPa;


        leafTurgorPsiCandidate = ...
            Sleaf.turgorPotential_MPa;


        leafApoplasticWaterCandidate = ...
            Sleaf.apoplasticWaterAmount_mol_m2;


        leafSymplasticWaterCandidate = ...
            Sleaf.symplasticWaterAmount_mol_m2;

    end


    if useLeafCapacitance

        if ...
                ~isfinite(leafStoragePrevious) ...
                || ~isfinite(dtSeconds) ...
                || dtSeconds <= 0

            state.rejectCode = "invalid_leaf_storage_previous";
            return

        end


        leafStorageRate = ...
            ( ...
                leafStorageAmountCandidate ...
                - leafStoragePrevious ...
            ) ...
            ./ dtSeconds;


    else

        % Initial timestep establishes the leaf storage state
        % without imposing a pre-simulation storage flux.

        leafStorageRate = ...
            0;

    end


else

    leafStorageAmountCandidate = ...
        NaN;


    leafCapacitanceCandidate = ...
        0;


    leafStorageBoundCandidate = ...
        false;


    leafStorageRate = ...
        0;

end

%% ========================================================================
% Leaf hydraulic network
% ========================================================================

flowTolerance = ...
    1e-8;


if leafNetCfg.enabled

    %% ------------------------------------------------------------
    % Symplastic storage rate
    %
    % The P-V model has fixed apoplastic water in the current
    % implementation, so:
    %
    %   dWleaf/dt = dWS/dt
    %
    % Nevertheless calculate dWS/dt explicitly from the symplastic
    % reservoir so the topology is clear.
    % ------------------------------------------------------------

    if ~isfinite(leafSymplasticWaterCandidate)

        state.rejectCode = "invalid_leaf_symplastic_state";
        return

    end


    if useLeafCapacitance

        if ...
                ~isfinite(leafSymplasticWaterPrevious) ...
                || ~isfinite(dtSeconds) ...
                || dtSeconds <= 0

            state.rejectCode = "invalid_leaf_symplastic_previous";
            return

        end


        symplasticStorageRate = ...
            ( ...
                leafSymplasticWaterCandidate ...
                - leafSymplasticWaterPrevious ...
            ) ...
            ./ dtSeconds;

    else

        symplasticStorageRate = ...
            0;

    end


    %% Difference should be essentially zero while apoplastic
    % water content is fixed.

    leafParallelStorageResidual = ...
        symplasticStorageRate ...
        - leafStorageRate;


    %% ------------------------------------------------------------
    % Solve parallel algebraic leaf network
    % ------------------------------------------------------------

    Lnet = ...
        leaf_parallel_hydraulic_network( ...
            psiLeaf, ...
            E, ...
            symplasticStorageRate, ...
            KXA, ...
            leafNetCfg.KapoplastSymplast, ...
            leafNetCfg.KsymplastEvaporation, ...
            leafNetCfg.KapoplastEvaporation);


    JleafIn = ...
        Lnet.JxylemToApoplast;


    psiX = ...
        Lnet.psiXylemRequired;

    if ...
            ~isfinite(JleafIn) ...
            || ~isfinite(psiX)

        state.rejectCode = "invalid_leaf_network_solution";
        return

    end


    if JleafIn < -flowTolerance

        if ~leafTraits.allowReverseXylemApoplastFlow

            state.rejectCode = "reverse_xylem_to_apoplast_flow";
            return

        end

    end


else

    %% ------------------------------------------------------------
    % Existing validated lumped-leaf formulation
    % ------------------------------------------------------------

    symplasticStorageRate = ...
        leafStorageRate;


    leafParallelStorageResidual = ...
        NaN;


    JleafIn = ...
        E ...
        + leafStorageRate;


    if JleafIn < -flowTolerance

        if ~leafTraits.allowReverseXylemApoplastFlow

            state.rejectCode = "reverse_xylem_to_apoplast_flow";
            return

        end

    end


    JleafIn = ...
        max(JleafIn,0);


    psiX = ...
        psiLeaf ...
        + JleafIn ./ Kleaf;


    Lnet = ...
        [];

end


if ~isfinite(psiX)

    state.rejectCode = "nonfinite_psi_xylem";
    return

end

if JleafIn < 0
    state.reverseLeafInternalFlow = JleafIn;
else
    state.reverseLeafInternalFlow = 0;
end

%% ========================================================================
% Xylem hydraulic storage
% ========================================================================

storageLumenBranchCodeCandidate = ...
    NaN;


storageLumenReleaseActiveCandidate = ...
    false;


storageLumenRefillActiveCandidate = ...
    false;


switch storageModel


    case "constantcapacitance"

        if useCapacitance

            if ...
                    ~isfinite(psiXPrevious) ...
                    || ~isfinite(dtSeconds) ...
                    || dtSeconds <= 0

                state.rejectCode = "invalid_xylem_pressure_previous";
                return

            end


            storageRate = ...
                Cx ...
                .* ( ...
                    psiX ...
                    - psiXPrevious ...
                ) ...
                ./ dtSeconds;

        else

            storageRate = ...
                0;

        end


    

        storageAmountCandidate = ...
            NaN;


        storageLumenCandidate = ...
            NaN;


        storageLivingCandidate = ...
            NaN;


        capacitanceEffectiveCandidate = ...
            Cx;


        capacitanceLumenCandidate = ...
            NaN;


        capacitanceLivingCandidate = ...
            NaN;


    case "twocomponentpiecewise"

        storageState = ...
            xylem_storage_pressure_curve( ...
                psiX, ...
                xylemTraits);


        storageAmountCandidate = ...
            storageState.total_mol_m2;


        storageLumenCandidate = ...
            storageState.lumen_mol_m2;


        storageLivingCandidate = ...
            storageState.living_mol_m2;


        capacitanceEffectiveCandidate = ...
            storageState.capacitanceTotal_mol_m2_MPa;


        capacitanceLumenCandidate = ...
            storageState.capacitanceLumen_mol_m2_MPa;


        capacitanceLivingCandidate = ...
            storageState.capacitanceLiving_mol_m2_MPa;


        if useCapacitance

            if ...
                    ~isfinite(storagePrevious) ...
                    || ~isfinite(dtSeconds) ...
                    || dtSeconds <= 0

                state.rejectCode = "invalid_xylem_storage_previous";
                return

            end


            storageRate = ...
                ( ...
                    storageAmountCandidate ...
                    - storagePrevious ...
                ) ...
                ./ dtSeconds;

        else

            storageRate = ...
                0;

        end


    case "twocomponenthysteretic"

        if useCapacitance

            hasPreviousStorage = ...
                   isfinite(psiXPrevious) ...
                && isfinite(storagePrevious) ...
                && isfinite(storageLumenPrevious);

        else

            hasPreviousStorage = ...
                false;

        end


        storageState = ...
            xylem_storage_hysteretic( ...
                psiX, ...
                psiXPrevious, ...
                storageLumenPrevious, ...
                xylemTraits, ...
                hasPreviousStorage);


        storageAmountCandidate = ...
            storageState.total_mol_m2;


        storageLumenCandidate = ...
            storageState.lumen_mol_m2;


        storageLivingCandidate = ...
            storageState.living_mol_m2;


        storageLumenBranchCodeCandidate = ...
            storageState.lumenBranchCode;


        storageLumenReleaseActiveCandidate = ...
            storageState.lumenReleaseActive;


        storageLumenRefillActiveCandidate = ...
            storageState.lumenRefillActive;


        capacitanceEffectiveCandidate = ...
            storageState.capacitanceTotal_mol_m2_MPa;


        capacitanceLumenCandidate = ...
            storageState.capacitanceLumen_mol_m2_MPa;


        capacitanceLivingCandidate = ...
            storageState.capacitanceLiving_mol_m2_MPa;


        if useCapacitance

            if ...
                    ~isfinite(storagePrevious) ...
                    || ~isfinite(dtSeconds) ...
                    || dtSeconds <= 0

                state.rejectCode = "invalid_xylem_hysteretic_previous";
                return

            end


            storageRate = ...
                ( ...
                    storageAmountCandidate ...
                    - storagePrevious ...
                ) ...
                ./ dtSeconds;

        else

            storageRate = ...
                0;

        end

end


%% ========================================================================
% Required root inflow from leaf demand and xylem storage
%
% Xylem node:
%
%   dSx/dt =
%       Jroot - JxylemToLeaf
%
% therefore:
%
%   Jroot =
%       JxylemToLeaf + dSx/dt
% ========================================================================

%% Diagnostic storage flux limiter

if isfield(xylemTraits,'diagnosticMaxStorageFlux')

    storageRate = ...
        max( ...
        min( ...
        storageRate, ...
        xylemTraits.diagnosticMaxStorageFlux), ...
        -xylemTraits.diagnosticMaxStorageFlux);

end

Jrequired = ...
    JleafIn ...
    + storageRate;

% A materially negative Jrequired represents net plant -> root/soil flow.
% Permit this only in the diagnostic reverse-root-flow case.
if Jrequired < -flowTolerance

    if ~leafTraits.allowReverseRootFlow

        state.rejectCode = "negative_required_root_flow";
        return

    end

elseif Jrequired < 0

    % Remove only tiny negative numerical noise.
    Jrequired = 0;

end

%% ========================================================================
% Nonlinear soil-root supply
% ========================================================================

%% ROOT DEMAND SOLVE TEST
rootBracketTolerance = 1e-10;

if isfield(rootPath,'rootDemandPsiMinMPa')
    psiRootMin = rootPath.rootDemandPsiMinMPa;
else
    psiRootMin = -5;
end

psiRootMax = 0;

rootDemandResidual = @(psiRoot) ...
    sum(soil_root_layer_flux( ...
        psiSoil,...
        psiRoot,...
        rootPath).flux,'omitnan') ...
    - Jrequired;


fRootMin = rootDemandResidual(psiRootMin);
fRootMax = rootDemandResidual(psiRootMax);

if ...
        ~isfinite(fRootMin) ...
        || ~isfinite(fRootMax)

    state.rejectCode = ...
        "nonfinite_root_demand_bracket";

    return

elseif abs(fRootMin) <= rootBracketTolerance

    psiRoot_demand = ...
        psiRootMin;

elseif abs(fRootMax) <= rootBracketTolerance

    psiRoot_demand = ...
        psiRootMax;

elseif sign(fRootMin) == sign(fRootMax)

    % Equal endpoint signs do not prove absence of an interior root because
    % the signed soil-root supply curve can be non-monotonic.
    %
    % Use the curve cached ONCE for this timestep to locate the preferred
    % least-negative interior sign-changing bracket. This avoids repeatedly
    % rescanning the same root-supply curve for every leaf-pressure
    % candidate.

    [interiorBracket,exactRoot,rootFound] = ...
        bracketRootDemandFromCacheLocal( ...
            rootDemandBracketCache, ...
            Jrequired, ...
            rootBracketTolerance, ...
            psiRootBranchReference);

    if ~rootFound

        state.rejectCode = ...
            "root_demand_unbracketed";

        return

    elseif isfinite(exactRoot)

        psiRoot_demand = ...
            exactRoot;

    else

        try

            psiRoot_demand = ...
                fzero( ...
                    rootDemandResidual, ...
                    interiorBracket);

        catch

            state.rejectCode = ...
                "root_demand_unbracketed";

            return

        end

    end

else

    psiRoot_demand = ...
        fzero( ...
            rootDemandResidual, ...
            [psiRootMin psiRootMax]);

end

% Root pressure is determined by root supply matching plant demand
psiRoot = psiRoot_demand;

layer = soil_root_layer_flux( ...
    psiSoil,...
    psiRoot,...
    rootPath);


if ...
        ~isfinite(psiRoot) ...
        || isempty(layer)

    state.rejectCode = "root_supply_failure";
    return

end

%% ========================================================================
% Dynamic stem conductance
% ========================================================================

Kstem = ...
    hydraulic_vulnerability( ...
        psiRoot, ...
        plantState.KstemMax, ...
        plantState.P50Stem, ...
        plantState.vulnerabilityShapeStem);


if ...
        ~isfinite(Kstem) ...
        || Kstem <= 0

    state.rejectCode = "invalid_stem_conductance";
    return

end


%% ========================================================================
% Actual stem inflow
% ========================================================================

stemGradient = psiRoot - psiX - gravityMPa;

Jstem = Kstem .* stemGradient;

Jbalance = Jstem - Jrequired;

if isempty(layer)
    % do nothing
else
    Jroot = sum(layer.flux,'omitnan');
end

state.rootDemandClosure = ...
    Jroot - Jrequired;

state.stemCapacityMargin = ...
    Jstem - Jrequired;

state.Jroot = Jroot;
state.Jstem = Jstem;
state.Jrequired = Jrequired;

residual = ...
    Jbalance;

%% ========================================================================
% Package candidate state
% ========================================================================

state.valid = ...
    true;

state.rejectCode = ...
    "valid";


state.E = ...
    E;


state.gs = ...
    gs;


state.psiRoot = ...
    psiRoot;


state.psiXylem = ...
    psiX;


state.Kstem = ...
    Kstem;


state.Kleaf = ...
    Kleaf;


%% Xylem storage

state.storageRate = ...
    storageRate;


state.storageAmount = ...
    storageAmountCandidate;


state.storageLumenAmount = ...
    storageLumenCandidate;


state.storageLivingAmount = ...
    storageLivingCandidate;


state.capacitanceEffective = ...
    capacitanceEffectiveCandidate;


state.capacitanceLumen = ...
    capacitanceLumenCandidate;


state.capacitanceLiving = ...
    capacitanceLivingCandidate;


state.storageLumenBranchCode = ...
    storageLumenBranchCodeCandidate;


state.storageLumenReleaseActive = ...
    storageLumenReleaseActiveCandidate;


state.storageLumenRefillActive = ...
    storageLumenRefillActiveCandidate;


%% Leaf storage

state.leafStorageRate = ...
    leafStorageRate;


state.leafStorageAmount = ...
    leafStorageAmountCandidate;


state.leafCapacitanceEffective = ...
    leafCapacitanceCandidate;


state.leafStorageBoundActive = ...
    leafStorageBoundCandidate;


%% Pressure-volume diagnostics

state.leafRelativeWaterContent = ...
    leafRWCCandidate;


state.leafSymplasticRelativeWaterContent = ...
    leafSymplasticRWCCandidate;


state.leafOsmoticPotential_MPa = ...
    leafOsmoticPsiCandidate;


state.leafTurgorPotential_MPa = ...
    leafTurgorPsiCandidate;


state.leafApoplasticWaterAmount_mol_m2 = ...
    leafApoplasticWaterCandidate;


state.leafSymplasticWaterAmount_mol_m2 = ...
    leafSymplasticWaterCandidate;

%% Parallel-network state

if leafNetCfg.enabled

    state.psiSymplast = ...
        psiLeaf;


    state.psiApoplast = ...
        Lnet.psiApoplast;


    state.psiEvaporationSite = ...
        Lnet.psiEvaporationSite;


    state.JxylemToApoplast = ...
        Lnet.JxylemToApoplast;


    state.JapoplastToSymplast = ...
        Lnet.JapoplastToSymplast;


    state.JsymplastToEvaporation = ...
        Lnet.JsymplastToEvaporation;


    state.JapoplastToEvaporation = ...
        Lnet.JapoplastToEvaporation;


    state.symplasticAdvectiveFraction = ...
        Lnet.symplasticAdvectiveFraction;


    state.apoplasticAdvectiveFraction = ...
        Lnet.apoplasticAdvectiveFraction;


    state.leafParallelApoplastClosure = ...
        Lnet.apoplastClosure;


    state.leafParallelSymplastClosure = ...
        Lnet.symplastClosure;


    state.leafParallelEvaporationClosure = ...
        Lnet.evaporationClosure;


    state.leafParallelWholeClosure = ...
        Lnet.wholeLeafClosure;


    state.leafParallelStorageRateResidual = ...
        leafParallelStorageResidual;


    state.KxylemApoplast = ...
        KXA;

end

state.JxylemToLeaf = ...
    JleafIn;


%% Stem

state.requiredStemFlow = ...
    Jrequired;


state.actualStemFlow = ...
    Jstem;

state.layer = ...
    layer;

end

%% =========================================================================
function S = leaf_storage_linear_state( ...
    psiLeaf, ...
    cfg)
%LEAF_STORAGE_LINEAR_STATE
%
% Linear development pressure-volume relation:
%
%   S =
%       Sref + C * (psiLeaf - psiRef)


Sref = ...
    cfg.referenceWaterContent_mol_m2;


psiRef = ...
    cfg.referencePsi_MPa;


C = ...
    cfg.capacitance_mol_m2_MPa;


rawAmount = ...
    Sref ...
    + C ...
    .* ( ...
        psiLeaf ...
        - psiRef);


amount = ...
    rawAmount;


boundActive = ...
    false;


if isfield(cfg,'minimumWaterContent_mol_m2')

    Wmin = ...
        cfg.minimumWaterContent_mol_m2;


    if amount < Wmin

        amount = ...
            Wmin;


        boundActive = ...
            true;

    end

end


if isfield(cfg,'maximumWaterContent_mol_m2')

    Wmax = ...
        cfg.maximumWaterContent_mol_m2;


    if amount > Wmax

        amount = ...
            Wmax;


        boundActive = ...
            true;

    end

end


if ...
        ~isfinite(amount) ...
        || amount <= 0

    error( ...
        'Leaf hydraulic storage became non-positive.');

end


if boundActive

    capacitanceEffective = ...
        0;

else

    capacitanceEffective = ...
        C;

end


S.amount_mol_m2 = ...
    amount;


S.rawAmount_mol_m2 = ...
    rawAmount;


S.capacitanceEffective_mol_m2_MPa = ...
    capacitanceEffective;


S.boundActive = ...
    boundActive;

end



%% ===============================================================
function value = getTimeValue(x, t)

if isscalar(x)

    value = x;

else

    value = x(t);

end

end


%% ===============================================================
function es = saturationVaporPressureKPa(TC)

es = ...
    0.61078 .* ...
    exp( ...
    (17.2694 .* TC) ./ ...
    (TC + 237.29));

end

function [psiR, layer] = solveRootPotential( ...
    psiS, ...
    pathModel, ...
    E, ...
    psiMin, ...
    allowReverseNetRootFlow)
%SOLVEROOTPOTENTIAL
%
% Solve:
%
%   sum_i J_i(psiRoot) = E
%
% using nonlinear matric-flux-potential soil/root fluxes.
%
% With hydraulic redistribution enabled, individual layer fluxes
% may be signed:
%
%   J_i > 0   soil -> root
%   J_i < 0   root -> soil
%
% while E remains the NET root-system inflow required by the plant.

psiS = ...
    psiS(:);


layer = [];


valid = ...
    isfinite(psiS);


%% ========================================================================
% Redistribution mode
% ========================================================================

allowHR = ...
    isfield(pathModel,'allowHydraulicRedistribution') ...
    && logical(pathModel.allowHydraulicRedistribution);


%% ========================================================================
% Input feasibility
%
% Negative E means net water movement from the plant/root system toward
% the soil. This requires both:
%
%   1. the diagnostic reverse-root-flow switch; and
%   2. a root pathway capable of signed hydraulic redistribution.
% ========================================================================

if ...
        ~any(valid) ...
        || ~isfinite(E)

    psiR = ...
        NaN;

    return

end


if ...
        E < 0 ...
        && ~(allowReverseNetRootFlow && allowHR)

    psiR = ...
        NaN;

    return

end


%% ========================================================================
% Historical zero-transpiration behavior
%
% Uptake-only mode:
%
%   put psiRoot at the wettest soil potential so every layer has zero
%   or positive driving potential and therefore zero realized uptake.
%
% Redistribution mode:
%
%   DO NOT special-case E = 0.
%
%   Instead solve:
%
%       sum_i J_i = 0
%
%   permitting simultaneous uptake and release.
% ========================================================================

if ...
        ~allowHR ...
        && E <= eps

    psiR = ...
        max(psiS(valid));


    layer = ...
        soil_root_layer_flux( ...
            psiS, ...
            psiR, ...
            pathModel);

    return

end


%% ========================================================================
% Root-potential bounds
% ========================================================================

psiUpper = max(psiSoil(valid)) + 0.05;

psiLower = psiMin;

psiUpper = max( ...
    psiUpper,...
    psiLower + 1e-6);


if psiLower >= psiUpper

    psiR = ...
        NaN;

    return

end


%% ========================================================================
% Supply at lower bound
% ========================================================================

Llower = ...
    soil_root_layer_flux( ...
        psiS, ...
        psiLower, ...
        pathModel);


supplyLower = ...
    sum( ...
        Llower.flux, ...
        'omitnan');


tolE = ...
    1e-10 ...
    .* max(1,abs(E));


if E > supplyLower + tolE

    psiR = ...
        NaN;

    return

end


%% ========================================================================
% Upper-bound residual
%
% In redistribution mode, psiUpper is the wettest soil potential.
%
% At this point:
%
%   wettest layer       J ~= 0
%   all drier layers    J <= 0
%
% so net signed supply should be <= 0.
% ========================================================================

Lupper = ...
    soil_root_layer_flux( ...
        psiS, ...
        psiUpper, ...
        pathModel);


resLower = ...
    supplyLower ...
    - E;


resUpper = ...
    sum( ...
        Lupper.flux, ...
        'omitnan') ...
    - E;


%% Exact-boundary solutions

if abs(resLower) <= tolE

    psiR = ...
        psiLower;


    layer = ...
        Llower;

    return

end


if abs(resUpper) <= tolE

    psiR = ...
        psiUpper;


    layer = ...
        Lupper;

    return

end


%% Require a genuine bracket

if ...
        ~isfinite(resLower) ...
        || ~isfinite(resUpper) ...
        || resLower .* resUpper > 0

    psiR = ...
        NaN;

    return

end


%% ========================================================================
% Supply residual
% ========================================================================

fun = @(p) ...
    rootSupplyResidual( ...
        p, ...
        psiS, ...
        pathModel, ...
        E);


%% ========================================================================
% Monotonic bracketed solve
% ========================================================================

options = ...
    optimset( ...
        'Display','off', ...
        'TolX',1e-10);


try

    psiR = ...
        fzero( ...
            fun, ...
            [psiLower psiUpper], ...
            options);

catch

    psiR = ...
        NaN;

    return

end


%% ========================================================================
% Final layer solution
% ========================================================================

layer = ...
    soil_root_layer_flux( ...
        psiS, ...
        psiR, ...
        pathModel);

end


%% ========================================================================
function r = rootSupplyResidual( ...
    psiRoot, ...
    psiS, ...
    pathModel, ...
    E)

L = soil_root_layer_flux( ...
    psiS, ...
    psiRoot, ...
    pathModel);

r = ...
    sum(L.flux,'omitnan') ...
    - E;

end

function [psiSolution, info] = ...
    findPhysicalHydraulicRoot( ...
    fun, psiMin, psiMax, psiPrevious)
%FINDPHYSICALHYDRAULICROOT
%
% Search the allowed hydraulic domain for roots of the coupled
% leaf-water-potential equation.
%
% The previous accepted state is searched first using progressively wider
% brackets. The original dense global/local scan is retained as a fallback.


%% ========================================================================
% Configuration
% ========================================================================

nGlobal = ...
    600;


% J2b outer-search candidate v2:
% Preserve the existing 600-point global fallback and the +/-0.5 MPa dense
% fallback. A progressive dense continuation search at +/-0.5, +/-1, and
% +/-2 MPa with ~0.001 MPa spacing is applied before the global fallback.
nLocal = ...
    1001;


localHalfWidth_MPa = ...
    0.5;


sampleRootTolerance = ...
    1e-10;


%% ========================================================================
% Validate hydraulic bounds
% ========================================================================

if ...
        ~isscalar(psiMin) ...
        || ~isscalar(psiMax) ...
        || ~isfinite(psiMin) ...
        || ~isfinite(psiMax) ...
        || psiMax <= psiMin

    error( ...
        'Invalid hydraulic search bounds.');

end


%% ========================================================================
% Global hydraulic search grid
% ========================================================================

psiGlobal = ...
    linspace( ...
        psiMin, ...
        psiMax, ...
        nGlobal);


usePrevious = ...
    isscalar(psiPrevious) ...
    && isfinite(psiPrevious);


%% ========================================================================
% Fast continuation-first adaptive bracket
% ========================================================================

if usePrevious

    psiPreviousSearch = ...
        min( ...
            max(psiPrevious,psiMin), ...
            psiMax);


    adaptiveHalfWidths_MPa = ...
        [ ...
        0.0025 ...
        0.005 ...
        0.01 ...
        0.02 ...
        0.05 ...
        0.10 ...
        0.25 ...
        0.50 ...
        1.0 ...
        2.0 ...
        4.0];


    psiEval = ...
        psiPreviousSearch;


    try

        fCenter = ...
            fun(psiPreviousSearch);

    catch

        fCenter = ...
            NaN;

    end


    fEval = ...
        fCenter;


    if ...
            isfinite(fCenter) ...
            && abs(fCenter) <= sampleRootTolerance

        psiSolution = ...
            psiPreviousSearch;


        info.rootFound = ...
            true;


        info.nRoots = ...
            1;


        info.allRoots = ...
            psiSolution;


        info.residual = ...
            fCenter;


        info.usedPreviousStateLocalGrid = ...
            true;


        info.previousStateSearchPsi = ...
            psiPreviousSearch;


        info.localHalfWidth_MPa = ...
            0;


        info.nGlobalGridPoints = ...
            nGlobal;


        info.nLocalGridPoints = ...
            1;


        info.nTotalSearchPoints = ...
            1;


        info.nFiniteSearchPoints = ...
            1;


        info.minimumAbsResidual = ...
            abs(fCenter);


        info.closestPsi = ...
            psiPreviousSearch;


        info.closestResidual = ...
            fCenter;


        info.usedAdaptiveBracket = ...
            true;


        info.adaptiveBracketHalfWidth_MPa = ...
            0;


        return

    end


    for kAdaptive = 1:numel(adaptiveHalfWidths_MPa)

        w = ...
            adaptiveHalfWidths_MPa(kAdaptive);


        psiLeft = ...
            max( ...
                psiMin, ...
                psiPreviousSearch - w);


        psiRight = ...
            min( ...
                psiMax, ...
                psiPreviousSearch + w);


        try

            fLeft = ...
                fun(psiLeft);

        catch

            fLeft = ...
                NaN;

        end


        try

            fRight = ...
                fun(psiRight);

        catch

            fRight = ...
                NaN;

        end


        psiEval = ...
            [psiEval psiLeft psiRight]; 


        fEval = ...
            [fEval fLeft fRight]; 


        brackets = ...
            zeros(0,2);


        if ...
                isfinite(fLeft) ...
                && abs(fLeft) <= sampleRootTolerance

            psiSolution = ...
                psiLeft;


            rootResidual = ...
                fLeft;


        elseif ...
                isfinite(fRight) ...
                && abs(fRight) <= sampleRootTolerance

            psiSolution = ...
                psiRight;


            rootResidual = ...
                fRight;


        else

            psiSolution = ...
                NaN;


            rootResidual = ...
                NaN;


            if ...
                    isfinite(fLeft) ...
                    && isfinite(fCenter) ...
                    && fLeft .* fCenter < 0

                brackets(end+1,:) = ...
                    [psiLeft psiPreviousSearch]; 

            end


            if ...
                    isfinite(fCenter) ...
                    && isfinite(fRight) ...
                    && fCenter .* fRight < 0

                brackets(end+1,:) = ...
                    [psiPreviousSearch psiRight]; 

            end


            if ...
                    isempty(brackets) ...
                    && isfinite(fLeft) ...
                    && isfinite(fRight) ...
                    && fLeft .* fRight < 0

                brackets(end+1,:) = ...
                    [psiLeft psiRight]; 

            end


            if ~isempty(brackets)

                candidateRoots = ...
                    NaN(1,size(brackets,1));


                candidateResiduals = ...
                    NaN(size(candidateRoots));


                for kb = 1:size(brackets,1)

                    try

                        rr = ...
                            fzero( ...
                                fun, ...
                                brackets(kb,:));


                        fr = ...
                            fun(rr);


                        if ...
                                isfinite(rr) ...
                                && isfinite(fr) ...
                                && rr >= psiMin ...
                                && rr <= psiMax

                            candidateRoots(kb) = ...
                                rr;


                            candidateResiduals(kb) = ...
                                fr;

                        end

                    catch

                        % Invalid bracket/root candidate.

                    end

                end


                validRoots = ...
                    isfinite(candidateRoots);


                if any(validRoots)

                    rootsTmp = ...
                        candidateRoots(validRoots);


                    [~,iiBest] = ...
                        min( ...
                            abs( ...
                                rootsTmp ...
                                - psiPreviousSearch));


                    psiSolution = ...
                        rootsTmp(iiBest);


                    resTmp = ...
                        candidateResiduals(validRoots);


                    rootResidual = ...
                        resTmp(iiBest);

                end

            end

        end


        if isfinite(psiSolution)

            finiteEval = ...
                isfinite(fEval);


            [minAbsTmp,iiTmp] = ...
                min( ...
                    abs( ...
                        fEval(finiteEval)));


            finitePsiTmp = ...
                psiEval(finiteEval);


            finiteFTmp = ...
                fEval(finiteEval);


            info.rootFound = ...
                true;


            info.nRoots = ...
                1;


            info.allRoots = ...
                psiSolution;


            info.residual = ...
                rootResidual;


            info.usedPreviousStateLocalGrid = ...
                true;


            info.previousStateSearchPsi = ...
                psiPreviousSearch;


            info.localHalfWidth_MPa = ...
                w;


            info.nGlobalGridPoints = ...
                nGlobal;


            info.nLocalGridPoints = ...
                numel(unique(psiEval));


            info.nTotalSearchPoints = ...
                numel(psiEval);


            info.nFiniteSearchPoints = ...
                sum(finiteEval);


            info.minimumAbsResidual = ...
                minAbsTmp;


            info.closestPsi = ...
                finitePsiTmp(iiTmp);


            info.closestResidual = ...
                finiteFTmp(iiTmp);


            info.usedAdaptiveBracket = ...
                true;


            info.adaptiveBracketHalfWidth_MPa = ...
                w;


            return

        end


        if ...
                psiLeft <= psiMin ...
                && psiRight >= psiMax

            break

        end

    end

end


%% ========================================================================
% Progressive dense continuation search
% ========================================================================
%
% The endpoint-only adaptive bracket above can jump over a narrow physical
% root. Before invoking the coarse global fallback, search progressively
% wider dense corridors around the previous accepted psiLeaf while retaining
% approximately 0.001 MPa spacing.
%
% Stop at the first corridor containing a sampled root or a finite
% sign-changing interval. This preserves the inexpensive +/-0.5 MPa search
% for ordinary states and only expands to +/-1 or +/-2 MPa when needed.

if usePrevious

    progressiveHalfWidths_MPa = ...
        [0.5 1.0 2.0];


    progressiveSpacing_MPa = ...
        0.001;


    for kProgressive = 1:numel(progressiveHalfWidths_MPa)

        wProgressive = ...
            progressiveHalfWidths_MPa(kProgressive);


        [psiProgressive, progressiveInfo] = ...
            findDenseLocalHydraulicRootLocal( ...
                fun, ...
                psiMin, ...
                psiMax, ...
                psiPreviousSearch, ...
                wProgressive, ...
                progressiveSpacing_MPa, ...
                sampleRootTolerance);


        if progressiveInfo.rootFound

            psiSolution = ...
                psiProgressive;


            info.rootFound = ...
                true;


            info.nRoots = ...
                progressiveInfo.nRoots;


            info.allRoots = ...
                progressiveInfo.allRoots;


            info.residual = ...
                progressiveInfo.residual;


            info.usedPreviousStateLocalGrid = ...
                true;


            info.previousStateSearchPsi = ...
                psiPreviousSearch;


            info.localHalfWidth_MPa = ...
                wProgressive;


            info.nGlobalGridPoints = ...
                nGlobal;


            info.nLocalGridPoints = ...
                progressiveInfo.nGridPoints;


            info.nTotalSearchPoints = ...
                progressiveInfo.nGridPoints;


            info.nFiniteSearchPoints = ...
                progressiveInfo.nFiniteGridPoints;


            info.minimumAbsResidual = ...
                progressiveInfo.minimumAbsResidual;


            info.closestPsi = ...
                progressiveInfo.closestPsi;


            info.closestResidual = ...
                progressiveInfo.closestResidual;


            info.usedAdaptiveBracket = ...
                false;


            info.adaptiveBracketHalfWidth_MPa = ...
                NaN;


            info.usedProgressiveDenseLocalGrid = ...
                true;


            info.progressiveDenseHalfWidth_MPa = ...
                wProgressive;


            return

        end

    end

end


%% ========================================================================
% Dense global/local fallback
% ========================================================================

if usePrevious

    psiPreviousSearch = ...
        min( ...
            max(psiPrevious,psiMin), ...
            psiMax);


    psiLocalMin = ...
        max( ...
            psiMin, ...
            psiPreviousSearch ...
            - localHalfWidth_MPa);


    psiLocalMax = ...
        min( ...
            psiMax, ...
            psiPreviousSearch ...
            + localHalfWidth_MPa);


    if psiLocalMax > psiLocalMin

        psiLocal = ...
            linspace( ...
                psiLocalMin, ...
                psiLocalMax, ...
                nLocal);

    else

        psiLocal = ...
            psiPreviousSearch;

    end


    psiGrid = ...
        unique( ...
            sort( ...
                [ ...
                psiGlobal(:); ...
                psiLocal(:); ...
                psiPreviousSearch]));

else

    psiPreviousSearch = ...
        NaN;


    psiLocal = ...
        [];


    psiGrid = ...
        psiGlobal(:);

end


psiGrid = ...
    psiGrid(:).';


nSearch = ...
    numel(psiGrid);


%% ========================================================================
% Evaluate hydraulic residual
% ========================================================================

F = ...
    NaN(size(psiGrid));


for i = 1:nSearch

    try

        F(i) = ...
            fun(psiGrid(i));

    catch

        F(i) = ...
            NaN;

    end

end


finiteMask = ...
    isfinite(F);


%% ========================================================================
% Minimum-residual diagnostic
% ========================================================================

if any(finiteMask)

    finiteIndices = ...
        find(finiteMask);


    [minimumAbsResidual,ii] = ...
        min( ...
            abs(F(finiteMask)));


    closestIndex = ...
        finiteIndices(ii);


    closestPsi = ...
        psiGrid(closestIndex);


    closestResidual = ...
        F(closestIndex);

else

    minimumAbsResidual = ...
        NaN;


    closestPsi = ...
        NaN;


    closestResidual = ...
        NaN;

end


%% ========================================================================
% Detect sampled roots
% ========================================================================

sampleRootMask = ...
    finiteMask ...
    & abs(F) <= sampleRootTolerance;


rootsFound = ...
    psiGrid(sampleRootMask);


%% ========================================================================
% Locate sign-changing intervals
% ========================================================================

validPair = ...
    isfinite(F(1:end-1)) ...
    & isfinite(F(2:end));


strictSignChange = ...
    validPair ...
    & ( ...
        F(1:end-1) ...
        .* F(2:end) ...
        < 0);


idx = ...
    find(strictSignChange);


%% ========================================================================
% Solve detected roots
% ========================================================================

for j = 1:numel(idx)

    bracket = [ ...
        psiGrid(idx(j)), ...
        psiGrid(idx(j)+1)];


    try

        r = ...
            fzero( ...
                fun, ...
                bracket);


        if ...
                isfinite(r) ...
                && r >= psiMin ...
                && r <= psiMax

            Fr = ...
                fun(r);


            if isfinite(Fr)

                rootsFound(end+1) = ...
                    r; 

            end

        end

    catch

        % Invalid hydraulic bracket.

    end

end


%% ========================================================================
% Remove duplicate roots
% ========================================================================

if ~isempty(rootsFound)

    rootsFound = ...
        sort( ...
            uniquetol( ...
                rootsFound, ...
                1e-8));

end


%% ========================================================================
% Select physical solution
% ========================================================================

if ~isempty(rootsFound)

    if usePrevious

        [~,best] = ...
            min( ...
                abs( ...
                    rootsFound ...
                    - psiPreviousSearch));

    else

        best = ...
            1;

    end


    psiSolution = ...
        rootsFound(best);


    info.rootFound = ...
        true;

else

    psiSolution = ...
        NaN;


    info.rootFound = ...
        false;

end


%% ========================================================================
% Root diagnostics
% ========================================================================

info.nRoots = ...
    numel(rootsFound);


info.allRoots = ...
    rootsFound;


if isfinite(psiSolution)

    try

        info.residual = ...
            fun(psiSolution);

    catch

        info.residual = ...
            NaN;

    end

else

    info.residual = ...
        NaN;

end


%% ========================================================================
% Search diagnostics
% ========================================================================

info.usedPreviousStateLocalGrid = ...
    usePrevious;


info.previousStateSearchPsi = ...
    psiPreviousSearch;


info.localHalfWidth_MPa = ...
    localHalfWidth_MPa;


info.nGlobalGridPoints = ...
    nGlobal;


info.nLocalGridPoints = ...
    numel(psiLocal);


info.nTotalSearchPoints = ...
    nSearch;


info.nFiniteSearchPoints = ...
    sum(finiteMask);


info.minimumAbsResidual = ...
    minimumAbsResidual;


info.closestPsi = ...
    closestPsi;


info.closestResidual = ...
    closestResidual;


info.usedAdaptiveBracket = ...
    false;


info.adaptiveBracketHalfWidth_MPa = ...
    NaN;

end

function S = leaf_storage_state(psiLeaf,cfg)
%LEAF_STORAGE_STATE
%
% Dispatch leaf hydraulic storage to the selected constitutive model.


if isfield(cfg,'model')

    modelName = ...
        lower(string(cfg.model));

else

    % Historical/default model.
    modelName = ...
        "linear";

end


switch modelName


    case "linear"

        S = ...
            leaf_storage_linear_state( ...
                psiLeaf, ...
                cfg);


    case { ...
            "pressurevolume", ...
            "pressure-volume", ...
            "pv"}

        S = ...
            leaf_storage_pressure_volume_state( ...
                psiLeaf, ...
                cfg);


    otherwise

        error( ...
            'Unknown leaf storage model: %s', ...
            modelName);

end

end


function KXA = resolveLeafXylemApoplastConductance( ...
    psiSymplast, ...
    plantState, ...
    leafNetCfg)
%RESOLVELEAFXYLEMAPOPLASTCONDUCTANCE


mode = ...
    lower(string( ...
        leafNetCfg.KxylemApoplastMode));


switch mode


    case "fixed"

        KXA = ...
            leafNetCfg.KxylemApoplast;


    case "legacyvulnerability"

        KXA = ...
            hydraulic_vulnerability( ...
                psiSymplast, ...
                plantState.KleafMax, ...
                plantState.P50Leaf, ...
                plantState.vulnerabilityShapeLeaf);


    otherwise

        error( ...
            'Unknown KxylemApoplastMode.');

end

end

%% ========================================================================


function C = buildRootDemandBracketCacheLocal( ...
    psiSoil, ...
    rootPath, ...
    psiMin, ...
    psiMax, ...
    nPoint)
%BUILDROOTDEMANDBRACKETCACHELOCAL
%
% Build one root-supply curve for a fixed soil state. This is shared by all
% leaf-pressure candidates at the timestep.

if nargin < 5 || isempty(nPoint)
    nPoint = 33;
end

psi = linspace(psiMax,psiMin,nPoint); % least-negative -> most-negative
supply = NaN(size(psi));

for k = 1:numel(psi)

    L = soil_root_layer_flux( ...
        psiSoil, ...
        psi(k), ...
        rootPath);

    supply(k) = ...
        sum(L.flux,'omitnan');

end

C = struct;
C.psiRoot_MPa = psi;
C.netRootSupply_mol_m2_s = supply;

end


function [bracket,exactRoot,found] = ...
    bracketRootDemandFromCacheLocal( ...
        C, ...
        Jrequired, ...
        tol, ...
        psiReference)
%BRACKETROOTDEMANDFROMCACHELOCAL
%
% Locate all resolved sign-changing root-demand brackets on the cached
% nonlinear root-supply curve, then select the branch closest to the
% previous accepted root pressure.
%
% If no finite continuation reference exists, fall back to the
% least-negative admissible crossing.

bracket = [NaN NaN];
exactRoot = NaN;
found = false;

psi = C.psiRoot_MPa(:).';
r = C.netRootSupply_mol_m2_s(:).' - Jrequired;

finiteMask = isfinite(psi) & isfinite(r);

iExactAll = find(finiteMask & abs(r) <= tol);

brackets = zeros(0,2);
bracketCenters = zeros(0,1);

for k = 1:(numel(psi)-1)

    if ...
            finiteMask(k) ...
            && finiteMask(k+1) ...
            && sign(r(k)) ~= sign(r(k+1))

        b = sort([psi(k) psi(k+1)]);

        brackets(end+1,:) = b; 
        bracketCenters(end+1,1) = mean(b); 

    end

end

if isempty(iExactAll) && isempty(brackets)
    return
end

if isfinite(psiReference)

    exactPsi = psi(iExactAll).';

    candidatePsi = [exactPsi; bracketCenters];
    candidateType = [ones(numel(exactPsi),1); ...
        2.*ones(numel(bracketCenters),1)];

    [~,j] = min(abs(candidatePsi-psiReference));

    if candidateType(j)==1

        exactRoot = candidatePsi(j);
        found = true;
        return

    else

        jBracket = j-numel(exactPsi);
        bracket = brackets(jBracket,:);
        found = true;
        return

    end

else

    candidatePsi = [];

    if ~isempty(iExactAll)
        candidatePsi = [candidatePsi; psi(iExactAll).']; 
    end

    if ~isempty(bracketCenters)
        candidatePsi = [candidatePsi; bracketCenters]; 
    end

    [~,j] = max(candidatePsi);

    nExact = numel(iExactAll);

    if j<=nExact
        exactRoot = candidatePsi(j);
    else
        bracket = brackets(j-nExact,:);
    end

    found = true;

end

end

function [psiSolution, info] = ...
    findDenseLocalHydraulicRootLocal( ...
    fun, psiMin, psiMax, psiPrevious, halfWidth_MPa, spacing_MPa, ...
    sampleRootTolerance)
%FINDDENSELOCALHYDRAULICROOTLOCAL
%
% Dense continuation-only hydraulic root search used by
% findPhysicalHydraulicRoot. The grid is centered on the previous accepted
% psiLeaf, clipped to the physical search domain, and resolved at
% approximately spacing_MPa.

psiLocalMin = ...
    max( ...
        psiMin, ...
        psiPrevious - halfWidth_MPa);


psiLocalMax = ...
    min( ...
        psiMax, ...
        psiPrevious + halfWidth_MPa);


if psiLocalMax > psiLocalMin

    nGrid = ...
        max( ...
            3, ...
            round( ...
                (psiLocalMax - psiLocalMin) ...
                / spacing_MPa) ...
            + 1);


    psiGrid = ...
        linspace( ...
            psiLocalMin, ...
            psiLocalMax, ...
            nGrid);

else

    psiGrid = ...
        psiPrevious;

end


F = ...
    NaN(size(psiGrid));


for i = 1:numel(psiGrid)

    try

        F(i) = ...
            fun(psiGrid(i));

    catch

        F(i) = ...
            NaN;

    end

end


finiteMask = ...
    isfinite(F);


if any(finiteMask)

    finiteIndices = ...
        find(finiteMask);


    [minimumAbsResidual,ii] = ...
        min(abs(F(finiteMask)));


    closestIndex = ...
        finiteIndices(ii);


    closestPsi = ...
        psiGrid(closestIndex);


    closestResidual = ...
        F(closestIndex);

else

    minimumAbsResidual = ...
        NaN;


    closestPsi = ...
        NaN;


    closestResidual = ...
        NaN;

end


sampleRootMask = ...
    finiteMask ...
    & abs(F) <= sampleRootTolerance;


rootsFound = ...
    psiGrid(sampleRootMask);


if numel(psiGrid) >= 2

    validPair = ...
        isfinite(F(1:end-1)) ...
        & isfinite(F(2:end));


    strictSignChange = ...
        validPair ...
        & ( ...
            F(1:end-1) ...
            .* F(2:end) ...
            < 0);


    idx = ...
        find(strictSignChange);

else

    idx = ...
        [];

end


for j = 1:numel(idx)

    bracket = ...
        [ ...
        psiGrid(idx(j)), ...
        psiGrid(idx(j)+1)];


    try

        r = ...
            fzero( ...
                fun, ...
                bracket);


        if ...
                isfinite(r) ...
                && r >= psiMin ...
                && r <= psiMax

            Fr = ...
                fun(r);


            if isfinite(Fr)

                rootsFound(end+1) = ...
                    r; 

            end

        end

    catch

        % Invalid local hydraulic bracket.

    end

end


%% ------------------------------------------------------------------------
% Refine finite/nonfinite admissibility edges
% ------------------------------------------------------------------------
edgeRefineSpacing_MPa = ...
    5e-5;


finiteTransition = ...
    find( ...
        finiteMask(1:end-1) ...
        ~= finiteMask(2:end));


for jEdge = 1:numel(finiteTransition)

    iEdge = ...
        finiteTransition(jEdge);


    psiEdgeA = ...
        psiGrid(iEdge);


    psiEdgeB = ...
        psiGrid(iEdge+1);


    edgeWidth = ...
        abs(psiEdgeB - psiEdgeA);


    nEdge = ...
        max( ...
            3, ...
            ceil(edgeWidth / edgeRefineSpacing_MPa) + 1);


    psiEdgeGrid = ...
        linspace(psiEdgeA,psiEdgeB,nEdge);


    Fedge = ...
        NaN(size(psiEdgeGrid));


    for qEdge = 1:numel(psiEdgeGrid)

        try
            Fedge(qEdge) = fun(psiEdgeGrid(qEdge));
        catch
            Fedge(qEdge) = NaN;
        end

    end


    finiteEdge = ...
        isfinite(Fedge);


    sampleEdgeRootMask = ...
        finiteEdge ...
        & abs(Fedge) <= sampleRootTolerance;


    if any(sampleEdgeRootMask)
        rootsFound = [rootsFound, psiEdgeGrid(sampleEdgeRootMask)]; 
    end


    if numel(psiEdgeGrid) >= 2

        validEdgePair = ...
            finiteEdge(1:end-1) ...
            & finiteEdge(2:end);


        strictEdgeSignChange = ...
            validEdgePair ...
            & (Fedge(1:end-1).*Fedge(2:end) < 0);


        idxEdge = ...
            find(strictEdgeSignChange);


        for qCross = 1:numel(idxEdge)

            kCross = idxEdge(qCross);

            bracketEdge = ...
                [psiEdgeGrid(kCross), psiEdgeGrid(kCross+1)];


            try

                rEdge = ...
                    fzero(fun,bracketEdge);


                if ...
                        isfinite(rEdge) ...
                        && rEdge >= psiMin ...
                        && rEdge <= psiMax

                    FrEdge = ...
                        fun(rEdge);


                    if isfinite(FrEdge)
                        rootsFound(end+1) = rEdge; 
                    end

                end

            catch

                % Invalid refined admissibility-edge bracket.

            end

        end

    end

end


if ~isempty(rootsFound)

    rootsFound = ...
        sort( ...
            uniquetol( ...
                rootsFound, ...
                1e-8));


    [~,best] = ...
        min( ...
            abs( ...
                rootsFound ...
                - psiPrevious));


    psiSolution = ...
        rootsFound(best);


    try

        residual = ...
            fun(psiSolution);

    catch

        residual = ...
            NaN;

    end


    rootFound = ...
        isfinite(residual);

else

    psiSolution = ...
        NaN;


    residual = ...
        NaN;


    rootFound = ...
        false;

end


info = ...
    struct;


info.rootFound = ...
    rootFound;


info.nRoots = ...
    numel(rootsFound);


info.allRoots = ...
    rootsFound;


info.residual = ...
    residual;


info.nGridPoints = ...
    numel(psiGrid);


info.nFiniteGridPoints = ...
    sum(finiteMask);


info.minimumAbsResidual = ...
    minimumAbsResidual;


info.closestPsi = ...
    closestPsi;


info.closestResidual = ...
    closestResidual;

end

