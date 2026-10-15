function out = intracellular_water_dynamic_exactD(timeDays,leafWater,cellTraits)
%INTRACELLULAR_WATER_DYNAMIC_EXACTD
%
% Physiologically parameterized exact-D intracellular water-isotope model.
%
% PRODUCTION ARCHITECTURE
%
%   validated bulk spatial symplast / mesophyll water
%                         |
%                         | gross tracer exchange + storage adjustment
%                         v
%                      cytosol
%                         |
%                         | gross tracer exchange + storage adjustment
%                         v
%                       plastid
%
% Optional:
%
%                      cytosol
%                         |
%                         v
%                   ER-associated water
%
%
% -------------------------------------------------------------------------
% PHYSIOLOGICAL PARAMETERIZATION
%
% Intracellular water amounts are fractions of the authoritative spatial
% symplast water amount:
%
%   W_Cy(t) = f_Cy * W_S(t)
%   W_P (t) = f_P  * W_S(t)
%   W_ER(t) = f_ER * W_S(t)       [optional]
%
% Gross nonfractionating tracer exchange uses interface-specific rates.
%
% Symplast -> cytosol remains an effective residence-time parameterization:
%
%   Q_SCy(t) = W_Cy(t) / tau_SCy
%
% Cytosol <-> plastid supports two backward-compatible modes:
%
%   legacyTimescale:
%       k_CyP = 1 / tau_CyP
%
%   membranePhysical:
%       k_CyP = Pm * (A_m / V_P)
%
% In both modes:
%
%   Q_CyP(t) = W_P(t) * k_CyP
%
% Optional cytosol <-> ER exchange remains:
%
%   Q_CyER(t) = W_ER(t) / tau_CyER
%
% Thus Q is a gross molecular exchange flux and may be much larger than
% net transpiration.
%
% -------------------------------------------------------------------------
% CHANGING WATER STORAGE
%
% All active intracellular pools scale with W_S(t). Therefore their water
% amounts may change between hydraulic nodes.
%
% For an interval:
%
%   s_Cy = dW_Cy/dt
%   s_P  = dW_P/dt
%   s_ER = dW_ER/dt
%
% The conservative storage-water network is:
%
%   S -> Cy : s_Cy + s_P + s_ER
%   Cy -> P : s_P
%   Cy -> ER: s_ER
%
% Positive flow is downstream; negative flow reverses donor direction.
% This construction guarantees the prescribed compartment water changes
% exactly while preserving donor D transport.
%
% -------------------------------------------------------------------------
% EXACT-D BOOKKEEPING
%
% Every active intracellular compartment carries:
%
%   ND_j = 2 * W_j * xD_j
%
% Gross exchange on an interface i <-> j is:
%
%   F_D,gross(i->j) = 2 * Q_ij * (xD_i - xD_j)
%
% Storage-water advection:
%
%   S_ij >= 0 : F_D,storage = 2*S_ij*xD_i
%   S_ij <  0 : F_D,storage = 2*S_ij*xD_j
%
%
% -------------------------------------------------------------------------
% REQUIRED INPUTS
%
% timeDays
%   Strictly increasing model times, days.
%
% leafWater
%   Isotope source fields:
%       ratioSymplastBulk
%       deltaDSymplastBulk
%
%   and symplast water amount:
%       symplasticWaterAmount_mol_m2
%
%   For synthetic tests, cellTraits.symplastWaterAmount_mol_m2 may be
%   supplied instead.
%
% cellTraits.waterFractionOfSymplast.cytosol
% cellTraits.waterFractionOfSymplast.plastid       if plastid enabled
%
% cellTraits.exchangeTimescale_s.symplastCytosol
%
% For plastid exchange, default legacy mode:
%   cellTraits.exchangeTimescale_s.cytosolPlastid
%
% Optional physical mode:
%   cellTraits.cytosolPlastidExchange.mode = 'membranePhysical'
%   cellTraits.cytosolPlastidExchange.membraneWaterPermeability_m_s
%   cellTraits.cytosolPlastidExchange.membraneAreaToPlastidWaterVolume_m_inv
%
% with:
%   k_CyP = Pm * (A_m/V_P)
%
% -------------------------------------------------------------------------
% RECOMMENDED FIRST-PASS PRODUCTION PRIORS
%
%   includePlastid = true
%   includeER = false
%
%   waterFractionOfSymplast.cytosol = 0.05
%   waterFractionOfSymplast.plastid = 0.003
%
%   exchangeTimescale_s.symplastCytosol = 0.5
%   exchangeTimescale_s.cytosolPlastid = 0.1
%
% Suggested exchange-timescale sensitivity:
%
%   [0.02 0.1 0.5 1 5 30 300] s
%
% -------------------------------------------------------------------------
% OPTIONAL INPUTS
%
% cellTraits.includePlastid         default true
% cellTraits.includeER              default false
%
% If ER enabled:
%   waterFractionOfSymplast.ER
%   exchangeTimescale_s.cytosolER
%
% cellTraits.sourcePool
%   'symplastBulk'       [default]
%   'apoplast'
%   'evaporationSite'
%   'symplastCell'
%
% cellTraits.sourceCellIndex        for sourcePool='symplastCell'
%
% cellTraits.sourceInterpolation
%   'linear'             [default]
%   'next'
%
% cellTraits.initialState
%   'sourceEqual'        [default]
%   'specified'
%
% For specified state:
%   initialDeltaD_permil.cytosol
%   initialDeltaD_permil.plastid
%   initialDeltaD_permil.ER
%
% or equivalent initialRatioNormalized fields.
%
% cellTraits.odeRelTol               default 1e-10
% cellTraits.odeAbsTol               default 1e-14
% cellTraits.odeMaxStepFraction      default 1
%
% cellTraits.R_VSMOW                 default 155.76e-6
%
% -------------------------------------------------------------------------


%% ========================================================================
% Time
% ========================================================================

t = timeDays(:);

if numel(t) < 2 || any(~isfinite(t)) || any(diff(t) <= 0)
    error('timeDays must be finite and strictly increasing.');
end

nTime = numel(t);
secondsPerDay = 86400;

%% ========================================================================
% VSMOW
% ========================================================================

if isfield(cellTraits,'R_VSMOW')
    Rstd = cellTraits.R_VSMOW;
else
    Rstd = 155.76e-6;
end

if ~isscalar(Rstd) || ~isfinite(Rstd) || Rstd <= 0
    error('cellTraits.R_VSMOW must be a positive finite scalar.');
end

%% ========================================================================
% Active intracellular compartments
% ========================================================================

if isfield(cellTraits,'includePlastid')
    includeP = logical(cellTraits.includePlastid);
else
    includeP = true;
end

if isfield(cellTraits,'includeER')
    includeER = logical(cellTraits.includeER);
else
    includeER = false;
end

if ~isscalar(includeP) || ~isscalar(includeER)
    error('includePlastid and includeER must be scalar logical values.');
end

% Explicitly reject the superseded chloroplast naming in the new
% production parameterization.
if isfield(cellTraits,'includeChloroplast')
    error([ ...
        'includeChloroplast has been superseded by includePlastid in the ', ...
        'wax-specific intracellular model.']);
end

names = {'cytosol'};

if includeP
    names{end+1} = 'plastid'; 
end

if includeER
    names{end+1} = 'ER'; 
end

nComp = numel(names);

iCy = find(strcmp(names,'cytosol'),1);
iP  = find(strcmp(names,'plastid'),1);
iER = find(strcmp(names,'ER'),1);

%% ========================================================================
% Reject superseded explicit W/Q parameterization
% ========================================================================

if isfield(cellTraits,'waterAmount_mol_m2') || ...
        isfield(cellTraits,'exchangeFlux_mol_m2_s')
    error([ ...
        'The production intracellular model now derives water amounts ', ...
        'from waterFractionOfSymplast and gross exchange from ', ...
        'exchangeTimescale_s. Explicit waterAmount_mol_m2 / ', ...
        'exchangeFlux_mol_m2_s inputs are retained only in archived ', ...
        'validation versions.']);
end

if isfield(cellTraits,'chloroplast') || ...
        (isfield(cellTraits,'waterFractionOfSymplast') && ...
         isfield(cellTraits.waterFractionOfSymplast,'chloroplast')) || ...
        (isfield(cellTraits,'exchangeTimescale_s') && ...
         isfield(cellTraits.exchangeTimescale_s,'cytosolChloroplast'))
    error('Use plastid / cytosolPlastid.');
end

%% ========================================================================
% Fractions of authoritative symplast water
% ========================================================================

if ~isfield(cellTraits,'waterFractionOfSymplast') || ...
        ~isstruct(cellTraits.waterFractionOfSymplast)
    error('cellTraits.waterFractionOfSymplast must be a structure.');
end

F = cellTraits.waterFractionOfSymplast;

fraction = zeros(nComp,1);

fraction(iCy) = requireFractionField( ...
    F,'cytosol','cellTraits.waterFractionOfSymplast');

if includeP
    fraction(iP) = requireFractionField( ...
        F,'plastid','cellTraits.waterFractionOfSymplast');
end

if includeER
    fraction(iER) = requireFractionField( ...
        F,'ER','cellTraits.waterFractionOfSymplast');
end

if sum(fraction) >= 1
    error([ ...
        'The sum of active intracellular water fractions must remain ', ...
        'below 1 because the compartments are subsets of symplast water.']);
end

%% ========================================================================
% Exchange parameterization
% ========================================================================

if ~isfield(cellTraits,'exchangeTimescale_s') || ...
        ~isstruct(cellTraits.exchangeTimescale_s)
    error('cellTraits.exchangeTimescale_s must be a structure.');
end

Tau = cellTraits.exchangeTimescale_s;

% Bulk symplast -> cytosol remains an effective access/mixing timescale.
tauSCy = requirePositiveOrInfField( ...
    Tau,'symplastCytosol','cellTraits.exchangeTimescale_s');

% Cytosol <-> plastid supports either the legacy target-pool timescale or
% the physical membrane relation k = Pm*(A/Vp).
if includeP

    if isfield(cellTraits,'cytosolPlastidExchange') && ...
            isstruct(cellTraits.cytosolPlastidExchange) && ...
            isfield(cellTraits.cytosolPlastidExchange,'mode')

        cyPMode = lower(char(cellTraits.cytosolPlastidExchange.mode));

    else

        cyPMode = 'legacytimescale';

    end

    switch cyPMode

        case 'legacytimescale'

            tauCyP = requirePositiveOrInfField( ...
                Tau,'cytosolPlastid','cellTraits.exchangeTimescale_s');

            if isinf(tauCyP)
                kCyP = 0;
            else
                kCyP = 1 ./ tauCyP;
            end

            PmCyP = NaN;
            areaToVolumeCyP = NaN;

        case 'membranephysical'

            X = cellTraits.cytosolPlastidExchange;

            PmCyP = requirePositiveField( ...
                X,'membraneWaterPermeability_m_s', ...
                'cellTraits.cytosolPlastidExchange');

            areaToVolumeCyP = requirePositiveField( ...
                X,'membraneAreaToPlastidWaterVolume_m_inv', ...
                'cellTraits.cytosolPlastidExchange');

            kCyP = PmCyP .* areaToVolumeCyP;

            if ~isscalar(kCyP) || ~isfinite(kCyP) || kCyP <= 0
                error('Derived cytosol-plastid k must be positive and finite.');
            end

            tauCyP = 1 ./ kCyP;

        otherwise

            error([ ...
                'cellTraits.cytosolPlastidExchange.mode must be ', ...
                '''legacyTimescale'' or ''membranePhysical''.']);

    end

else

    cyPMode = 'disabled';
    tauCyP = inf;
    kCyP = 0;
    PmCyP = NaN;
    areaToVolumeCyP = NaN;

end

if includeER
    tauCyER = requirePositiveOrInfField( ...
        Tau,'cytosolER','cellTraits.exchangeTimescale_s');
else
    tauCyER = inf;
end

%% ========================================================================
% Authoritative symplast water amount
% ========================================================================

WS = selectSymplastWaterAmount(leafWater,cellTraits,nTime);

if any(~isfinite(WS)) || any(WS <= 0)
    error('Symplast water amount must remain positive and finite.');
end

Wnode = fraction .* WS.';

%% ========================================================================
% Leaf-water isotope source
% ========================================================================

if isfield(cellTraits,'sourcePool')
    sourcePool = lower(char(cellTraits.sourcePool));
else
    sourcePool = 'symplastbulk';
end

switch sourcePool

    case 'symplastbulk'

        [rSource,inputField] = selectRatioSeries( ...
            leafWater, ...
            {'ratioSymplastBulk','ratioSymplast'}, ...
            {'deltaDSymplastBulk','deltaDSymplast'}, ...
            nTime);

        sourceDescription = ...
            'bulk spatial symplast / mesophyll water';

    case 'apoplast'

        [rSource,inputField] = selectRatioSeries( ...
            leafWater, ...
            {'ratioApoplast'}, ...
            {'deltaDApoplast'}, ...
            nTime);

        sourceDescription = 'apoplastic water';

    case 'evaporationsite'

        [rSource,inputField] = selectRatioSeries( ...
            leafWater, ...
            {'ratioEvaporationSite'}, ...
            {'deltaDEvaporationSite'}, ...
            nTime);

        sourceDescription = 'evaporation-site water';

    case 'symplastcell'

        if ~isfield(cellTraits,'sourceCellIndex')
            error('sourceCellIndex is required for sourcePool = symplastCell.');
        end

        idx = cellTraits.sourceCellIndex;

        if ~isscalar(idx) || ~isfinite(idx) || idx < 1 || idx ~= round(idx)
            error('sourceCellIndex must be a positive integer.');
        end

        if isfield(leafWater,'ratioSymplastCell')

            rSource = extractCellSeries( ...
                leafWater.ratioSymplastCell,idx,nTime);

            inputField = 'ratioSymplastCell';

        elseif isfield(leafWater,'deltaDSymplastCell')

            dSource = extractCellSeries( ...
                leafWater.deltaDSymplastCell,idx,nTime);

            rSource = 1 + dSource./1000;

            inputField = 'deltaDSymplastCell';

        else
            error([ ...
                'symplastCell source requires ratioSymplastCell or ', ...
                'deltaDSymplastCell.']);
        end

        sourceDescription = sprintf('spatial symplast cell %d',idx);

    otherwise

        error('Unknown sourcePool: %s',sourcePool);

end

rSource = rSource(:);

if numel(rSource) ~= nTime || any(~isfinite(rSource)) || any(rSource <= 0)
    error('Selected leaf-water isotope source is invalid.');
end

RSource = Rstd .* rSource;
xSource = RSource ./ (1 + RSource);
deltaSource = 1000 .* (rSource - 1);

%% ========================================================================
% Source interpolation
% ========================================================================

if isfield(cellTraits,'sourceInterpolation')
    sourceInterpolation = lower(char(cellTraits.sourceInterpolation));
else
    sourceInterpolation = 'linear';
end

if ~ismember(sourceInterpolation,{'linear','next'})
    error('sourceInterpolation must be ''linear'' or ''next''.');
end

%% ========================================================================
% Initial isotope state
% ========================================================================

if isfield(cellTraits,'initialState')
    initialState = lower(char(cellTraits.initialState));
else
    initialState = 'sourceequal';
end

r0 = zeros(nComp,1);

switch initialState

    case 'sourceequal'

        r0(:) = rSource(1);

    case 'specified'

        if isfield(cellTraits,'initialRatioNormalized')

            S0 = cellTraits.initialRatioNormalized;

            for j = 1:nComp
                r0(j) = requirePositiveField( ...
                    S0,names{j}, ...
                    'cellTraits.initialRatioNormalized');
            end

        elseif isfield(cellTraits,'initialDeltaD_permil')

            S0 = cellTraits.initialDeltaD_permil;

            for j = 1:nComp

                d0 = requireFiniteField( ...
                    S0,names{j}, ...
                    'cellTraits.initialDeltaD_permil');

                r0(j) = 1 + d0./1000;

            end

        else

            error([ ...
                'specified initialState requires initialRatioNormalized ', ...
                'or initialDeltaD_permil.']);

        end

    otherwise

        error('Unknown initialState: %s',initialState);

end

R0 = Rstd .* r0;
x0 = R0 ./ (1 + R0);

ND = zeros(nComp,nTime);
ND(:,1) = 2 .* Wnode(:,1) .* x0;

%% ========================================================================
% Solver configuration
% ========================================================================

if isfield(cellTraits,'odeRelTol')
    odeRelTol = cellTraits.odeRelTol;
else
    odeRelTol = 1e-10;
end

if isfield(cellTraits,'odeAbsTol')
    odeAbsTol = cellTraits.odeAbsTol;
else
    odeAbsTol = 1e-14;
end

if isfield(cellTraits,'odeMaxStepFraction')
    odeMaxStepFraction = cellTraits.odeMaxStepFraction;
else
    odeMaxStepFraction = 1;
end

if ~isscalar(odeRelTol) || ~isfinite(odeRelTol) || odeRelTol <= 0
    error('odeRelTol must be a positive finite scalar.');
end

if ~isscalar(odeAbsTol) || ~isfinite(odeAbsTol) || odeAbsTol <= 0
    error('odeAbsTol must be a positive finite scalar.');
end

if ~isscalar(odeMaxStepFraction) || ...
        ~isfinite(odeMaxStepFraction) || ...
        odeMaxStepFraction <= 0 || odeMaxStepFraction > 1
    error('odeMaxStepFraction must lie in (0,1].');
end

%% ========================================================================
% Cumulative exact-D transfer diagnostics
%
% C(1) = total source -> cytosol D transfer
% C(2) = total cytosol -> plastid D transfer
% C(3) = total cytosol -> ER D transfer
%
% Separate gross-exchange cumulative transfers are also tracked so the
% storage-driven component can be diagnosed independently.
% ========================================================================

Ctotal = zeros(3,nTime);
Cgross = zeros(3,nTime);

%% ========================================================================
% Time-varying derived water / exchange diagnostics
% ========================================================================

Qnode = zeros(3,nTime);

Qnode(1,:) = waterOverTau(Wnode(iCy,:),tauSCy);

if includeP
    Qnode(2,:) = Wnode(iP,:) .* kCyP;
end

if includeER
    Qnode(3,:) = waterOverTau(Wnode(iER,:),tauCyER);
end

storageWaterFlux = zeros(3,nTime);

%% ========================================================================
% Interval-wise exact-D integration
%
% Hydraulic/storage states W_S are linearly interpolated between accepted
% nodes. Intracellular fractions therefore vary linearly too.
%
% The gross exchange timescales remain fixed physical parameters while Q
% changes with target-pool water amount.
% ========================================================================

for k = 1:(nTime-1)

    dt = (t(k+1)-t(k)) .* secondsPerDay;

    W0 = Wnode(:,k);
    W1 = Wnode(:,k+1);
    dWdt = (W1-W0) ./ dt;

    % Signed storage water flows, positive in downstream direction.
    S_CyP = 0;
    S_CyER = 0;

    if includeP
        S_CyP = dWdt(iP);
    end

    if includeER
        S_CyER = dWdt(iER);
    end

    S_SCy = sum(dWdt);

    storageWaterFlux(:,k+1) = [S_SCy;S_CyP;S_CyER];

    switch sourceInterpolation

        case 'linear'
            xs0 = xSource(k);
            xs1 = xSource(k+1);

        case 'next'
            xs0 = xSource(k+1);
            xs1 = xSource(k+1);

    end

    % State:
    %   1:nComp       exact-D inventories
    %   next 3        cumulative total transfers
    %   next 3        cumulative gross transfers

    y0 = zeros(nComp+6,1);

    y0(1:nComp) = ND(:,k);

    y0(nComp+(1:3)) = Ctotal(:,k);
    y0(nComp+3+(1:3)) = Cgross(:,k);

    maxStep = max(dt .* odeMaxStepFraction,eps);

    opts = odeset( ...
        'RelTol',odeRelTol, ...
        'AbsTol',odeAbsTol, ...
        'MaxStep',maxStep);

    [~,Y] = ode15s( ...
        @(tt,yy) rhsInterval( ...
            tt,yy,dt, ...
            W0,W1,dWdt, ...
            xs0,xs1, ...
            tauSCy,kCyP,tauCyER, ...
            includeP,includeER, ...
            iCy,iP,iER, ...
            S_SCy,S_CyP,S_CyER), ...
        [0 dt], ...
        y0, ...
        opts);

    y1 = Y(end,:).';

    ND(:,k+1) = y1(1:nComp);

    Ctotal(:,k+1) = y1(nComp+(1:3));
    Cgross(:,k+1) = y1(nComp+3+(1:3));

end

%% ========================================================================
% Reconstruct exact isotope states
% ========================================================================

x = ND ./ (2 .* Wnode);

if any(~isfinite(x(:))) || any(x(:) <= 0) || any(x(:) >= 1)
    error('Intracellular D atom fraction left the physical interval (0,1).');
end

R = x ./ (1-x);
r = R ./ Rstd;
deltaD = 1000 .* (r-1);

%% ========================================================================
% Exact-D closure
% ========================================================================

deltaND = ND - ND(:,1);

expected = zeros(size(ND));

expected(iCy,:) = ...
    Ctotal(1,:) ...
    - Ctotal(2,:) ...
    - Ctotal(3,:);

if includeP
    expected(iP,:) = Ctotal(2,:);
end

if includeER
    expected(iER,:) = Ctotal(3,:);
end

compClosure = deltaND - expected;

totalND = sum(ND,1);

wholeClosure = ...
    (totalND-totalND(1)) ...
    - Ctotal(1,:);

Cstorage = Ctotal-Cgross;

%% ========================================================================
% Water-balance diagnostics
% ========================================================================

waterExpected = zeros(nComp,nTime);

waterExpected(:,1) = Wnode(:,1);

for k = 1:(nTime-1)

    dt = (t(k+1)-t(k)) .* secondsPerDay;

    S_SCy = storageWaterFlux(1,k+1);
    S_CyP = storageWaterFlux(2,k+1);
    S_CyER = storageWaterFlux(3,k+1);

    waterExpected(iCy,k+1) = ...
        waterExpected(iCy,k) ...
        + dt .* (S_SCy-S_CyP-S_CyER);

    if includeP
        waterExpected(iP,k+1) = ...
            waterExpected(iP,k) ...
            + dt .* S_CyP;
    end

    if includeER
        waterExpected(iER,k+1) = ...
            waterExpected(iER,k) ...
            + dt .* S_CyER;
    end

end

waterClosure = Wnode-waterExpected;

%% ========================================================================
% Outputs
% ========================================================================

out = struct;

out.timeDays = t.';

out.sourcePool = sourcePool;
out.sourcePoolDescription = sourceDescription;
out.sourceInputField = inputField;
out.sourceInterpolation = sourceInterpolation;

out.sourceRatioWater = rSource.';
out.sourceDeltaDWater = deltaSource.';
out.sourceDAtomFraction = xSource.';

out.symplasticWaterAmount_mol_m2 = WS.';

out.activeCompartments = names;
out.includePlastid = includeP;
out.includeER = includeER;

out.waterFractionOfSymplast = struct;
out.waterAmount_mol_m2 = struct;
out.exchangeTimescale_s = struct;
out.exchangeFlux_mol_m2_s = struct;

out.ratio = struct;
out.deltaD = struct;
out.DAtomFraction = struct;
out.DInventory_molD_m2 = struct;
out.compartmentDClosureResidual_molD_m2 = struct;
out.compartmentWaterClosureResidual_mol_m2 = struct;
out.apparentMemoryOffset_permil = struct;

for j = 1:nComp

    nm = names{j};

    out.waterFractionOfSymplast.(nm) = fraction(j);
    out.waterAmount_mol_m2.(nm) = Wnode(j,:);

    out.ratio.(nm) = r(j,:);
    out.deltaD.(nm) = deltaD(j,:);
    out.DAtomFraction.(nm) = x(j,:);
    out.DInventory_molD_m2.(nm) = ND(j,:);

    out.compartmentDClosureResidual_molD_m2.(nm) = ...
        compClosure(j,:);

    out.compartmentWaterClosureResidual_mol_m2.(nm) = ...
        waterClosure(j,:);

    out.apparentMemoryOffset_permil.(nm) = ...
        1000 .* (r(j,:)./rSource.' - 1);

end

out.exchangeTimescale_s.symplastCytosol = tauSCy;
out.exchangeTimescale_s.cytosolPlastid = tauCyP;
out.exchangeTimescale_s.cytosolER = tauCyER;

out.exchangeRate_s_inv = struct;
if isinf(tauSCy)
    out.exchangeRate_s_inv.symplastCytosol = 0;
else
    out.exchangeRate_s_inv.symplastCytosol = 1 ./ tauSCy;
end
out.exchangeRate_s_inv.cytosolPlastid = kCyP;
if isinf(tauCyER)
    out.exchangeRate_s_inv.cytosolER = 0;
else
    out.exchangeRate_s_inv.cytosolER = 1 ./ tauCyER;
end

out.cytosolPlastidExchange = struct;
out.cytosolPlastidExchange.mode = cyPMode;
out.cytosolPlastidExchange.k_s_inv = kCyP;
out.cytosolPlastidExchange.derivedTimescale_s = tauCyP;
out.cytosolPlastidExchange.membraneWaterPermeability_m_s = PmCyP;
out.cytosolPlastidExchange.membraneAreaToPlastidWaterVolume_m_inv = ...
    areaToVolumeCyP;

out.exchangeFlux_mol_m2_s.symplastCytosol = Qnode(1,:);
out.exchangeFlux_mol_m2_s.cytosolPlastid = Qnode(2,:);
out.exchangeFlux_mol_m2_s.cytosolER = Qnode(3,:);

out.storageWaterFlux_mol_m2_s = struct;
out.storageWaterFlux_mol_m2_s.symplastCytosol = storageWaterFlux(1,:);
out.storageWaterFlux_mol_m2_s.cytosolPlastid = storageWaterFlux(2,:);
out.storageWaterFlux_mol_m2_s.cytosolER = storageWaterFlux(3,:);

out.cumulativeDTransfer_molD_m2 = struct;
out.cumulativeDTransfer_molD_m2.symplastCytosol = Ctotal(1,:);
out.cumulativeDTransfer_molD_m2.cytosolPlastid = Ctotal(2,:);
out.cumulativeDTransfer_molD_m2.cytosolER = Ctotal(3,:);

out.cumulativeGrossExchangeDTransfer_molD_m2 = struct;
out.cumulativeGrossExchangeDTransfer_molD_m2.symplastCytosol = Cgross(1,:);
out.cumulativeGrossExchangeDTransfer_molD_m2.cytosolPlastid = Cgross(2,:);
out.cumulativeGrossExchangeDTransfer_molD_m2.cytosolER = Cgross(3,:);

out.cumulativeStorageDTransfer_molD_m2 = struct;
out.cumulativeStorageDTransfer_molD_m2.symplastCytosol = Cstorage(1,:);
out.cumulativeStorageDTransfer_molD_m2.cytosolPlastid = Cstorage(2,:);
out.cumulativeStorageDTransfer_molD_m2.cytosolER = Cstorage(3,:);

out.totalDInventory_molD_m2 = totalND;
out.wholeSystemDClosureResidual_molD_m2 = wholeClosure;

out.initialState = initialState;
out.R_VSMOW = Rstd;

out.intrinsicFractionationApplied = false;
out.alphaIntrinsic = 1;
out.epsilonIntrinsic_permil = 0;

out.biochemicalHydrogenSourcesApplied = false;
out.metabolicWaterApplied = false;

out.integrationMethod = 'intervalWiseStiffExactD';
out.parameterization = ...
    'symplastWaterFractions_plus_mixed_exchange_parameterization';

out.isotopeBookkeeping = 'exactDAtomInventory';
out.isHydraulicStorage = false;

out.modelRole = ...
    ['nested intracellular isotope-tracer network within the validated ', ...
     'leaf symplast; water amounts scale with authoritative symplast ', ...
     'storage but do not feed back on hydraulics'];

out.physicalExchangeUpgrade = ...
    ['cytosol-plastid exchange may be derived from membrane Pm*A/Vp; ', ...
     'symplast-cytosol remains an effective access/mixing timescale'];

end

%% =========================================================================
function dydt = rhsInterval( ...
    tt,y,dt, ...
    W0,W1,dWdt, ...
    xs0,xs1, ...
    tauSCy,kCyP,tauCyER, ...
    includeP,includeER, ...
    iCy,iP,iER, ...
    S_SCy,S_CyP,S_CyER)

nComp = numel(W0);

f = min(max(tt./dt,0),1);

W = W0 + f.*(W1-W0);

x = y(1:nComp) ./ (2.*W);

xSource = xs0 + f.*(xs1-xs0);

QSCy = waterOverTau(W(iCy),tauSCy);

if includeP
    QCyP = W(iP) .* kCyP;
else
    QCyP = 0;
end

if includeER
    QCyER = waterOverTau(W(iER),tauCyER);
else
    QCyER = 0;
end

% Gross nonfractionating tracer exchange.
Fg_SCy = 2 .* QSCy .* (xSource-x(iCy));

if includeP
    Fg_CyP = 2 .* QCyP .* (x(iCy)-x(iP));
else
    Fg_CyP = 0;
end

if includeER
    Fg_CyER = 2 .* QCyER .* (x(iCy)-x(iER));
else
    Fg_CyER = 0;
end

% Donor-aware D carried by changing compartment storage.
Fs_SCy = signedDAdvection(S_SCy,xSource,x(iCy));

if includeP
    Fs_CyP = signedDAdvection(S_CyP,x(iCy),x(iP));
else
    Fs_CyP = 0;
end

if includeER
    Fs_CyER = signedDAdvection(S_CyER,x(iCy),x(iER));
else
    Fs_CyER = 0;
end

F_SCy = Fg_SCy + Fs_SCy;
F_CyP = Fg_CyP + Fs_CyP;
F_CyER = Fg_CyER + Fs_CyER;

dND = zeros(nComp,1);

dND(iCy) = ...
    F_SCy ...
    - F_CyP ...
    - F_CyER;

if includeP
    dND(iP) = F_CyP;
end

if includeER
    dND(iER) = F_CyER;
end

dydt = zeros(nComp+6,1);

dydt(1:nComp) = dND;

% Total cumulative transfers.
dydt(nComp+1) = F_SCy;
dydt(nComp+2) = F_CyP;
dydt(nComp+3) = F_CyER;

% Gross-exchange-only cumulative transfers.
dydt(nComp+4) = Fg_SCy;
dydt(nComp+5) = Fg_CyP;
dydt(nComp+6) = Fg_CyER;

% dWdt is supplied explicitly for transparency and validation; the
% conservative signed storage flows above are constructed from it.
if any(~isfinite(dWdt))
    error('Nonfinite intracellular storage derivative encountered.');
end

end

%% =========================================================================
function FD = signedDAdvection(S,xUp,xDown)

if S >= 0
    FD = 2 .* S .* xUp;
else
    FD = 2 .* S .* xDown;
end

end

%% =========================================================================
function Q = waterOverTau(W,tau)

if isinf(tau)
    Q = zeros(size(W));
else
    Q = W ./ tau;
end

end

%% =========================================================================
function WS = selectSymplastWaterAmount(leafWater,cellTraits,nTime)

if isfield(leafWater,'symplasticWaterAmount_mol_m2')

    WS = leafWater.symplasticWaterAmount_mol_m2(:);

elseif isfield(leafWater,'symplasticCellWaterAmount_mol_m2')

    M = leafWater.symplasticCellWaterAmount_mol_m2;

    if isvector(M)

        if numel(M) ~= nTime
            error('symplasticCellWaterAmount_mol_m2 has incompatible dimensions.');
        end

        WS = M(:);

    elseif size(M,1) == nTime

        WS = sum(M,2);

    elseif size(M,2) == nTime

        WS = sum(M,1).';

    else

        error('symplasticCellWaterAmount_mol_m2 has incompatible dimensions.');

    end

elseif isfield(cellTraits,'symplastWaterAmount_mol_m2')

    v = cellTraits.symplastWaterAmount_mol_m2;

    if isscalar(v)
        WS = v .* ones(nTime,1);
    else
        WS = v(:);
    end

else

    error([ ...
        'Authoritative symplast water amount was not found. Supply ', ...
        'leafWater.symplasticWaterAmount_mol_m2 or, for synthetic tests, ', ...
        'cellTraits.symplastWaterAmount_mol_m2.']);

end

if numel(WS) ~= nTime
    error('Symplast water amount must contain one value per timestep.');
end

end

%% =========================================================================
function value = requireFractionField(S,name,parent)

if ~isstruct(S) || ~isfield(S,name)
    error('%s.%s is required.',parent,name);
end

value = S.(name);

if ~isscalar(value) || ~isfinite(value) || value <= 0 || value >= 1
    error('%s.%s must lie strictly between 0 and 1.',parent,name);
end

end

%% =========================================================================
function value = requirePositiveOrInfField(S,name,parent)

if ~isstruct(S) || ~isfield(S,name)
    error('%s.%s is required.',parent,name);
end

value = S.(name);

if ~isscalar(value) || isnan(value) || value <= 0
    error('%s.%s must be positive or Inf.',parent,name);
end

end

%% =========================================================================
function value = requirePositiveField(S,name,parent)

if ~isstruct(S) || ~isfield(S,name)
    error('%s.%s is required.',parent,name);
end

value = S.(name);

if ~isscalar(value) || ~isfinite(value) || value <= 0
    error('%s.%s must be a positive finite scalar.',parent,name);
end

end

%% =========================================================================
function value = requireFiniteField(S,name,parent)

if ~isstruct(S) || ~isfield(S,name)
    error('%s.%s is required.',parent,name);
end

value = S.(name);

if ~isscalar(value) || ~isfinite(value)
    error('%s.%s must be a finite scalar.',parent,name);
end

end

%% =========================================================================
function [r,inputField] = selectRatioSeries( ...
    S,ratioFields,deltaFields,nTime)

for i = 1:numel(ratioFields)

    f = ratioFields{i};

    if isfield(S,f)

        r = S.(f)(:);
        inputField = f;

        if numel(r) ~= nTime
            error('leafWater.%s must contain one value per timestep.',f);
        end

        return

    end

end

for i = 1:numel(deltaFields)

    f = deltaFields{i};

    if isfield(S,f)

        d = S.(f)(:);
        inputField = f;

        if numel(d) ~= nTime
            error('leafWater.%s must contain one value per timestep.',f);
        end

        r = 1 + d./1000;

        return

    end

end

error('Required leaf-water isotope source fields were not found.');

end

%% =========================================================================
function r = extractCellSeries(M,idx,nTime)

sz = size(M);

if isvector(M)

    if idx ~= 1 || numel(M) ~= nTime
        error('Spatial-cell isotope array has incompatible dimensions.');
    end

    r = M(:);
    return

end

if sz(1) == nTime

    if idx > sz(2)
        error('sourceCellIndex exceeds number of spatial cells.');
    end

    r = M(:,idx);

elseif sz(2) == nTime

    if idx > sz(1)
        error('sourceCellIndex exceeds number of spatial cells.');
    end

    r = M(idx,:).';

else

    error('Spatial-cell isotope array must have one dimension equal to nTime.');

end

end
