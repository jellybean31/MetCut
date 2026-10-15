function out = leaf_hydraulics_dynamic_spatial_symplast( ...
    timeDays, ...
    xylemWater, ...
    hydro, ...
    met, ...
    leafTraits)
%LEAF_HYDRAULICS_DYNAMIC_SPATIAL_SYMPLAST
%
% Dynamic exact-D leaf-water isotope model with a finite apoplastic
% reservoir and a spatially resolved symplastic / mesophyll pathway.
%
% Architecture:
%
%   xylem -> apoplast -> [S1 <-> S2 <-> ... <-> SN] -> evaporation site
%                  \----------------------------------> evaporation site
%                                                       <-> atmosphere
%
% The hydraulic water states are authoritative:
%
%   WA(t) = hydro.leafApoplasticWaterAmount_mol_m2
%   WS(t) = hydro.leafSymplasticWaterAmount_mol_m2
%
%   Wi(t) = WS(t)/N,  sum_i Wi(t) = WS(t).
%
% For hydraulic interval [t(k),t(k+1)], ending-node fluxes are held
% constant (FD convention):
%
%   Jx  = JxylemToApoplast(k+1)
%   JAS = JapoplastToSymplast(k+1)
%   JSE = JsymplastToEvaporation(k+1)
%   JAE = JapoplastToEvaporation(k+1)
%
% while WA and WS vary linearly between accepted hydraulic endpoints.
%
% The symplastic face-water flux is distributed so every equal-volume
% spatial cell has the exact hydraulic storage tendency:
%
%   Jface(i+1/2) = JAS - (i/N)*(JAS-JSE), i=0...N.
%
% Signed Scharfetter-Gummel advection-diffusion is used on internal faces
% and at the evaporation-site boundary. This naturally supports
% two-sided rehydration (A -> S(z) <- E) when JAS>0 and JSE<0.
%
% The evaporation site is massless. At every RHS evaluation its exact-D
% composition is solved simultaneously with the spatial right-boundary
% flux and fractionating atmospheric exchange. The direct apoplastic
% branch is retained. 
%
% Exact D-atom inventories are the prognostic isotope states.
%
% Spatial configuration (optional):
%   leafTraits.spatialSymplast.nCells              default 50
%   leafTraits.spatialSymplast.odeRelTol            default 1e-10
%   leafTraits.spatialSymplast.odeAbsTol            default 1e-15 mol D m^-2
%   leafTraits.spatialSymplast.odeMaxStepFraction   default 0.25
%   leafTraits.spatialSymplast.initialRatioCell     scalar or N values
%   leafTraits.spatialSymplast.initialDeltaDCell_permil scalar or N values
%   leafTraits.spatialSymplast.waterClosureTolerance_mol_m2_s default 1e-10

%% ========================================================================
% Time / constants
% ========================================================================

t = timeDays(:);
nTime = numel(t);

if nTime < 2 || any(~isfinite(t)) || any(diff(t) <= 0)
    error('timeDays must be a finite increasing vector.');
end

secondsPerDay = 86400;
R_VSMOW = 155.76e-6;

%% ========================================================================
% Current integration target requires the validated parallel leaf network
% ========================================================================

parallelHydraulics = ...
    isfield(hydro,'leafParallelHydraulicsEnabled') ...
    && logical(hydro.leafParallelHydraulicsEnabled);

if ~parallelHydraulics
    error([ ...
        'leaf_hydraulics_dynamic_spatial_symplast currently requires ' ...
        'hydro.leafParallelHydraulicsEnabled = true.']);
end

requiredHydro = { ...
    'leafApoplasticWaterAmount_mol_m2', ...
    'leafSymplasticWaterAmount_mol_m2', ...
    'leafHydraulicStorageRate_mol_m2_s', ...
    'E', ...
    'gs', ...
    'JxylemToApoplast', ...
    'JapoplastToSymplast', ...
    'JsymplastToEvaporation', ...
    'JapoplastToEvaporation'};

for k = 1:numel(requiredHydro)
    if ~isfield(hydro,requiredHydro{k})
        error('hydro.%s is required.',requiredHydro{k});
    end
end

WA = hydro.leafApoplasticWaterAmount_mol_m2(:);
WS = hydro.leafSymplasticWaterAmount_mol_m2(:);
dWLdtHydro = hydro.leafHydraulicStorageRate_mol_m2_s(:);
E = hydro.E(:);
gs = hydro.gs(:);
Jx = hydro.JxylemToApoplast(:);
JAS = hydro.JapoplastToSymplast(:);
JSE = hydro.JsymplastToEvaporation(:);
JAE = hydro.JapoplastToEvaporation(:);

allHydro = {WA,WS,dWLdtHydro,E,gs,Jx,JAS,JSE,JAE};
for k = 1:numel(allHydro)
    if numel(allHydro{k}) ~= nTime
        error('Hydraulic leaf arrays must match timeDays.');
    end
end

if any(~isfinite(WA) | WA <= 0) || any(~isfinite(WS) | WS <= 0)
    error('Apoplastic and symplastic water amounts must remain positive.');
end

if any(~isfinite([E;gs;Jx;JAS;JSE;JAE]))
    error('Leaf hydraulic flux/conductance arrays must remain finite.');
end

if any(gs < 0)
    error('hydro.gs must remain nonnegative.');
end

if isfield(leafTraits,'throughflowTolerance_mol_m2_s')

    throughflowTolerance = ...
        leafTraits.throughflowTolerance_mol_m2_s;

else

    throughflowTolerance = ...
        1e-8;

end

%% ========================================================================
% Xylem isotope boundary
% ========================================================================

if ~isfield(xylemWater,'ratioWater')
    error('xylemWater.ratioWater is required.');
end

rX = xylemWater.ratioWater(:);
if numel(rX) ~= nTime
    error('xylemWater.ratioWater must match timeDays.');
end

if any(~isfinite(rX) | rX <= 0)
    error('xylemWater.ratioWater must be finite and positive.');
end

if isfield(xylemWater,'Jwater_mol_m2_s')
    JxBoundary = xylemWater.Jwater_mol_m2_s(:);
    if numel(JxBoundary) ~= nTime
        error('xylemWater.Jwater_mol_m2_s must match timeDays.');
    end
    boundaryResidual = JxBoundary - Jx;
    if max(abs(boundaryResidual),[],'omitnan') > throughflowTolerance

        error( ...
            ['Xylem isotope boundary does not match hydraulic ' ...
            'JxylemToApoplast within the configured throughflow tolerance.']);

    end
else
    boundaryResidual = zeros(nTime,1);
end

%% ========================================================================
% Hydraulic closure / interval convention
% ========================================================================

dWAdtAmount = zeros(nTime,1);
dWSdtAmount = zeros(nTime,1);

for k = 2:nTime

    dt_s = ...
        (t(k)-t(k-1))*secondsPerDay;

    dWAdtAmount(k) = ...
        (WA(k)-WA(k-1))/dt_s;

    dWSdtAmount(k) = ...
        (WS(k)-WS(k-1))/dt_s;

end



dWSdtHydro = ...
    JAS - JSE;


apoplastWaterClosure = ...
    Jx - JAS - JAE - dWAdtAmount;

symplastWaterClosure = ...
    JAS - JSE - dWSdtHydro;

wholeLeafWaterClosure = ...
    Jx - E - dWLdtHydro;

pathwayFlowResidual = ...
    E - JAE - JSE;

hydraulicPathwayResidual = ...
    Jx - JAS - JAE;

intervalApoplastForcingResidual = NaN(nTime,1);
intervalSymplastForcingResidual = NaN(nTime,1);
intervalApoplastForcingResidual(2:end) = apoplastWaterClosure(2:end);
intervalSymplastForcingResidual(2:end) = symplastWaterClosure(2:end);

maxHydraulicResidual = max(abs([ ...
    apoplastWaterClosure(2:end); ...
    symplastWaterClosure(2:end); ...
    wholeLeafWaterClosure(2:end); ...
    pathwayFlowResidual(2:end); ...
    hydraulicPathwayResidual(2:end)]),[],'omitnan');

if maxHydraulicResidual > 1e-9
    error('Hydraulic pathway inputs are internally inconsistent.');
end

%% ========================================================================
% Spatial configuration
% ========================================================================
if isfield(leafTraits,'spatialSymplast')
    Scfg = leafTraits.spatialSymplast;
else
    Scfg = struct;
end

nCells = getCfgScalar(Scfg,'nCells',50);
if ~isfinite(nCells) || nCells < 1 || nCells ~= round(nCells)
    error('leafTraits.spatialSymplast.nCells must be a positive integer.');
end

odeRelTol = getCfgScalar(Scfg,'odeRelTol',1e-10);
odeAbsTol = getCfgScalar(Scfg,'odeAbsTol',1e-15);
odeMaxStepFraction = getCfgScalar(Scfg,'odeMaxStepFraction',0.25);
waterClosureTolerance = getCfgScalar( ...
    Scfg,'waterClosureTolerance_mol_m2_s',1e-10);

if odeRelTol <= 0 || odeAbsTol <= 0
    error('Spatial symplast ODE tolerances must be positive.');
end
if odeMaxStepFraction <= 0 || odeMaxStepFraction > 1
    error('spatialSymplast.odeMaxStepFraction must lie in (0,1).');
end
if waterClosureTolerance < 0
    error('spatialSymplast.waterClosureTolerance_mol_m2_s must be nonnegative.');
end

if any(abs(intervalSymplastForcingResidual(2:end)) > waterClosureTolerance)
    error([ ...
        'Symplast hydraulic forcing is inconsistent with the prescribed ' ...
        'water-storage trajectory under the ending-node convention.']);
end

%% ========================================================================
% Anatomy / liquid-water physical properties
% ========================================================================

if ~isfield(leafTraits,'anatomy')
    error('leafTraits.anatomy is required.');
end

Aleaf = leaf_anatomy_to_L(leafTraits.anatomy);

if abs(Aleaf.flowFraction - 1) > 1e-12
    error([ ...
        'The spatial symplast integration currently requires ' ...
        'leafTraits.anatomy.flowFraction = 1.']);
end

pathLength = Aleaf.pathLength_m;
aEffBase = Aleaf.effectiveFlowAreaFraction;

if ~isscalar(pathLength) || ~isfinite(pathLength) || pathLength <= 0
    error('Anatomy-derived symplastic path length must be positive.');
end
if ~isscalar(aEffBase) || ~isfinite(aEffBase) || aEffBase <= 0
    error('Anatomy-derived effective flow area must be positive.');
end

TleafC = expandTimeLocal(met.TleafC,nTime);

if isfield(leafTraits,'Dliquid_m2_s')
    DHDO = expandTimeLocal(leafTraits.Dliquid_m2_s,nTime);
else
    DHDO = hdoLiquidDiffusivityLocal(TleafC);
end

if any(~isfinite(DHDO) | DHDO <= 0)
    error('Liquid-water HDO diffusivity must be positive.');
end

if isfield(leafTraits,'waterMolarDensity_mol_m3')
    Cwater = expandTimeLocal(leafTraits.waterMolarDensity_mol_m3,nTime);
else
    Cwater = 55.5e3.*ones(nTime,1);
end

if any(~isfinite(Cwater) | Cwater <= 0)
    error('Liquid-water molar concentration must be positive.');
end

%% ========================================================================
% Symplastic diffusive participation / sensitivity scaling
% ========================================================================

PXcfg = struct;
if isfield(leafTraits,'outsideXylemPathways')
    PXcfg = leafTraits.outsideXylemPathways;
end

if isfield(PXcfg,'symplasticDiffusiveFraction')
    fSymDiff = expandTimeLocal(PXcfg.symplasticDiffusiveFraction,nTime);
else
    fSymDiff = ones(nTime,1);
end

if any(~isfinite(fSymDiff) | fSymDiff <= 0 | fSymDiff > 1)
    error([ ...
        'The spatial symplast model currently requires ' ...
        'symplasticDiffusiveFraction in (0,1).']);
end

if isfield(leafTraits,'isotopeDiffusionConductanceScale')
    KdiffScale = leafTraits.isotopeDiffusionConductanceScale;
else
    KdiffScale = 1;
end

if ~isscalar(KdiffScale) || ~isfinite(KdiffScale) || KdiffScale <= 0
    error('isotopeDiffusionConductanceScale must be positive.');
end

% Whole-path diffusive transport coefficient in mol H2O m^-2 s^-1 * m.
Dtilde = Cwater .* DHDO .* aEffBase .* fSymDiff .* KdiffScale;
KdiffSymWholePath = Dtilde ./ pathLength;

dz = pathLength/nCells;
zCenter_m = ((1:nCells).'-0.5).*dz;
zFace_m = (0:nCells).'.*dz;

%% ========================================================================
% Direct apoplast -> evaporation-site isotope diffusion
% ========================================================================

KdiffApo = zeros(nTime,1);

%% ========================================================================
% Meteorological / fractionation forcing
% ========================================================================

TairC = expandTimeLocal(met.TairC,nTime);
RH = expandTimeLocal(met.RH,nTime);
deltaDVapor = expandTimeLocal(met.deltaDVapor,nTime);

if isfield(met,'pressureKPa')
    pressure = expandTimeLocal(met.pressureKPa,nTime);
else
    pressure = 101.325.*ones(nTime,1);
end

if isfield(met,'gb_mol_m2_s')
    gb = expandTimeLocal(met.gb_mol_m2_s,nTime);
else
    gb = ones(nTime,1);
end

if any(gb <= 0)
    error('Boundary-layer conductance must be positive.');
end

rV = 1 + deltaDVapor./1000;

% Established liquid-vapor equilibrium fractionation.
TK = TleafC+273.15;
lnAlpha1000 = ...
    24.844e6./TK.^2 ...
    - 76.248e3./TK ...
    + 52.612;
alphaEq = exp(lnAlpha1000./1000);

% Resistance-weighted kinetic fractionation in zero-safe conductance form.
if ~isfield(leafTraits,'epsilonKStomatal') ...
        || ~isfield(leafTraits,'epsilonKBoundary')
    error('leafTraits.epsilonKStomatal and epsilonKBoundary are required.');
end

epsS = leafTraits.epsilonKStomatal;
epsB = leafTraits.epsilonKBoundary;
epsilonK = (epsS.*gb+epsB.*gs)./(gb+gs);
epsilonKRatio = epsilonK./1000;
alphaK = 1+epsilonKRatio;

if isfield(leafTraits,'disableIsotopeFractionation') ...
        && logical(leafTraits.disableIsotopeFractionation)
    alphaEq(:) = 1;
    alphaK(:) = 1;
end

rVin = rV./alphaK;

% Gross atmospheric water exchange, constrained to the authoritative net E.
ei = saturationVaporPressureKPaLocal(TleafC);
ea = RH .* saturationVaporPressureKPaLocal(TairC);
gw = gs.*gb./(gs+gb);
wi = ei./pressure;
wa = ea./pressure;
JoutNominal = gw.*wi;
JinNominal = gw.*wa;
JnetNominal = JoutNominal-JinNominal;
Jout = JoutNominal;
Jin = JinNominal;

for k = 1:nTime
    if k == 1 || mod(k,12) == 0 || k == nTime-1
        fprintf( ...
            'Leaf isotope interval %3d / %3d | %.3f -> %.3f d\n', ...
            k,nTime-1,t(k),t(k+1));
    end
    if abs(JnetNominal(k)) > 1e-14
        scale = E(k)/JnetNominal(k);
        Jout(k) = JoutNominal(k)*scale;
        Jin(k) = JinNominal(k)*scale;
    elseif abs(E(k)) > 1e-12
        error('Atmospheric exchange cannot reproduce nonzero E at timestep %d.',k);
    end
end

atmosphericWaterClosure = Jout-Jin-E;

%% ========================================================================
% Initial exact-D isotope state
% ========================================================================

rA0 = rX(1);
xA0 = ratioNormalizedToDAtomFraction(rA0,R_VSMOW);
NDA0 = 2*WA(1)*xA0;

if isfield(Scfg,'initialRatioCell')
    rCell0 = expandScalarOrCellVector( ...
        Scfg.initialRatioCell,nCells,'spatialSymplast.initialRatioCell');
elseif isfield(Scfg,'initialDeltaDCell_permil')
    dCell0 = expandScalarOrCellVector( ...
        Scfg.initialDeltaDCell_permil,nCells, ...
        'spatialSymplast.initialDeltaDCell_permil');
    rCell0 = 1 + dCell0./1000;
else
    rCell0 = rA0.*ones(nCells,1);
end

if any(~isfinite(rCell0) | rCell0 <= 0)
    error('Initial spatial symplast isotope ratios must be positive.');
end

xCell0 = ratioNormalizedToDAtomFraction(rCell0,R_VSMOW);
Wcell0 = WS(1)/nCells;
NDCell0 = 2.*Wcell0.*xCell0;

%% ========================================================================
% ODE state layout
%
% 1                           apoplast D inventory
% 2 : 1+N                     spatial symplast D inventories
% 2+N                         cumulative xylem D input
% 3+N                         cumulative atmospheric D input
% 4+N                         cumulative atmospheric D output
% 5+N                         cumulative direct A<->E D transfer
% 6+N : 6+N+N                 cumulative symplast face D transfers (N+1)
% ========================================================================

idxA = 1;
idxCell = 2:(1+nCells);
idxCumX = 2+nCells;
idxCumAtmIn = 3+nCells;
idxCumAtmOut = 4+nCells;
idxCumAE = 5+nCells;
idxCumFace = (6+nCells):(6+2*nCells);

nState = 6 + 2*nCells;
y0 = zeros(nState,1);
y0(idxA) = NDA0;
y0(idxCell) = NDCell0;

YSol = NaN(nTime,nState);
YSol(1,:) = y0.';
yCurrent = y0;

%% ========================================================================
% Interval-wise exact-D integration
% ========================================================================

for k = 1:(nTime-1)
    iForce = k+1;

    interval = struct;
    interval.t0 = t(k);
    interval.t1 = t(k+1);
    interval.WA0 = WA(k);
    interval.WA1 = WA(k+1);
    interval.WS0 = WS(k);
    interval.WS1 = WS(k+1);

    interval.rX = rX(iForce);
    interval.rVin = rVin(iForce);
    interval.Jx = JxBoundary(iForce);
    interval.JAS = JAS(iForce);
    interval.JSE = JSE(iForce);
    interval.JAE = JAE(iForce);
    interval.Jin = Jin(iForce);
    interval.Jout = Jout(iForce);
    interval.alphaEq = alphaEq(iForce);
    interval.alphaK = alphaK(iForce);
    interval.KdiffApo = KdiffApo(iForce);

    interval.Dtilde = Dtilde(iForce);
    interval.Ginternal = interval.Dtilde/dz;
    interval.Gright = 2*interval.Ginternal;
    interval.Jface = buildFaceWaterFlux(interval.JAS,interval.JSE,nCells);

    pInternal = interval.Jface(2:nCells)./interval.Ginternal;
    pRight = interval.Jface(end)./interval.Gright;
    interval.BminusInternal = bernoulliFunction(-pInternal);
    interval.BplusInternal = bernoulliFunction(pInternal);
    interval.BminusRight = bernoulliFunction(-pRight);
    interval.BplusRight = bernoulliFunction(pRight);

    interval.nCells = nCells;
    interval.idxA = idxA;
    interval.idxCell = idxCell;
    interval.idxCumX = idxCumX;
    interval.idxCumAtmIn = idxCumAtmIn;
    interval.idxCumAtmOut = idxCumAtmOut;
    interval.idxCumAE = idxCumAE;
    interval.idxCumFace = idxCumFace;

    interval.Rstd = R_VSMOW;
    interval.secondsPerDay = secondsPerDay;

    dtIntervalDays = interval.t1-interval.t0;

    odeAbsTolVector = odeAbsTol.*ones(nState,1);

    cumulativeIdx = [ ...
        idxCumX, ...
        idxCumAtmIn, ...
        idxCumAtmOut, ...
        idxCumAE, ...
        idxCumFace];

    odeAbsTolVector(cumulativeIdx) = 3e-12;

    odeOptions = odeset( ...
        'RelTol',odeRelTol, ...
        'AbsTol',odeAbsTolVector, ...
        'MaxStep',max(dtIntervalDays*odeMaxStepFraction,eps));

    [~,Yinterval] = ode15s( ...
        @(tt,yy) spatialLeafRHS(tt,yy,interval), ...
        [interval.t0 interval.t1], ...
        yCurrent, ...
        odeOptions);

    if isempty(Yinterval)
        error('Spatial leaf isotope integration failed on interval %d.',k);
    end

    yCurrent = Yinterval(end,:).';
    YSol(k+1,:) = yCurrent.';
end

%% ========================================================================
% Reconstruct finite isotope states
% ========================================================================

NDA = YSol(:,idxA);
NDCell = YSol(:,idxCell).';                 % N x time
Wcell = repmat((WS.'/nCells),nCells,1);

xA = NDA./(2.*WA);
xCell = NDCell./(2.*Wcell);

if any(~isfinite(xA) | xA <= 0 | xA >= 1) ...
        || any(~isfinite(xCell(:)) | xCell(:) <= 0 | xCell(:) >= 1)
    error('Integrated spatial leaf isotope state became nonphysical.');
end

rA = dAtomFractionToNormalizedRatio(xA,R_VSMOW);
rCell = dAtomFractionToNormalizedRatio(xCell,R_VSMOW);
deltaDA = 1000.*(rA-1);
deltaDCell = 1000.*(rCell-1);

NDsym = sum(NDCell,1).';
xSymBulk = NDsym./(2.*WS);
rSymBulk = dAtomFractionToNormalizedRatio(xSymBulk,R_VSMOW);
deltaDSymBulk = 1000.*(rSymBulk-1);

%% ========================================================================
% Node diagnostics: massless evaporation site and spatial face fluxes
% ========================================================================

rEvap = NaN(nTime,1);
xEvap = NaN(nTime,1);
faceWaterFlux = NaN(nCells+1,nTime);
faceDFlux = NaN(nCells+1,nTime);
siteDClosure = NaN(nTime,1);
directAEDFlux = NaN(nTime,1);

for k = 1:nTime
    Ginternal_k = Dtilde(k)/dz;
    Gright_k = 2*Ginternal_k;
    Jface_k = buildFaceWaterFlux(JAS(k),JSE(k),nCells);
    pRight_k = Jface_k(end)/Gright_k;
    BminusRight_k = bernoulliFunction(-pRight_k);
    BplusRight_k = bernoulliFunction(pRight_k);

    [rEvap(k),xEvap(k)] = solveSpatialEvaporationSiteExact( ...
        rA(k), ...
        xCell(end,k), ...
        rVin(k), ...
        JAE(k), ...
        Jin(k), ...
        Jout(k), ...
        KdiffApo(k), ...
        Gright_k, ...
        BminusRight_k, ...
        BplusRight_k, ...
        alphaEq(k), ...
        alphaK(k), ...
        R_VSMOW);

    faceWaterFlux(:,k) = Jface_k;
    faceFlux_xD_k = spatialFaceFluxesSignedPrecomputed( ...
        xCell(:,k), ...
        xA(k), ...
        xEvap(k), ...
        Jface_k, ...
        Ginternal_k, ...
        Gright_k, ...
        bernoulliFunction(-Jface_k(2:nCells)./Ginternal_k), ...
        bernoulliFunction(Jface_k(2:nCells)./Ginternal_k), ...
        BminusRight_k, ...
        BplusRight_k);

    faceDFlux(:,k) = 2.*faceFlux_xD_k;

    JAin = max(JAE(k),0);
    JAout = max(-JAE(k),0);
    FDAE = 2.*( ...
        (JAin+KdiffApo(k)).*xA(k) ...
        - (JAout+KdiffApo(k)).*xEvap(k));
    directAEDFlux(k) = FDAE;

    xVin_k = ratioNormalizedToDAtomFraction(rVin(k),R_VSMOW);
    rVout_k = rEvap(k)/(alphaEq(k)*alphaK(k));
    xVout_k = ratioNormalizedToDAtomFraction(rVout_k,R_VSMOW);
    FDatmIn_k = 2*Jin(k)*xVin_k;
    FDatmOut_k = 2*Jout(k)*xVout_k;

    siteDClosure(k) = ...
        FDAE ...
        + faceDFlux(end,k) ...
        + FDatmIn_k ...
        - FDatmOut_k;
end

deltaDEvap = 1000.*(rEvap-1);

%% ========================================================================
% Exact-D closure diagnostics from integrated cumulative transfers
% ========================================================================

cumX = YSol(:,idxCumX);
cumAtmIn = YSol(:,idxCumAtmIn);
cumAtmOut = YSol(:,idxCumAtmOut);
cumAE = YSol(:,idxCumAE);
cumFace = YSol(:,idxCumFace).';             % (N+1) x time

NDtotal = NDA + NDsym;
NDtotal0 = NDtotal(1);
DInventoryClosureResidual = ...
    (NDtotal-NDtotal0) ...
    - (cumX+cumAtmIn-cumAtmOut);

apoplastDClosureResidual = ...
    (NDA-NDA(1)) ...
    - (cumX-cumFace(1,:).'-cumAE);

cellDClosureResidual = NaN(nCells,nTime);
for i = 1:nCells
    cellDClosureResidual(i,:) = ...
        (NDCell(i,:)-NDCell(i,1)) ...
        - (cumFace(i,:)-cumFace(i+1,:));
end

symplastDClosureResidual = ...
    (NDsym-NDsym(1)) ...
    - (cumFace(1,:).'-cumFace(end,:).');

%% ========================================================================
% Spatial/path diagnostics
% ========================================================================

PeSymBoundary = JSE./KdiffSymWholePath;
phiClassical = NaN(nTime,1);
smallPe = abs(PeSymBoundary) < 1e-10;
phiClassical(smallPe) = 1;
phiClassical(~smallPe) = ...
    (1-exp(-PeSymBoundary(~smallPe)))./PeSymBoundary(~smallPe);

phiSpatialObserved = NaN(nTime,1);
den = xEvap-xA;
validPhi = abs(den) > 1e-16;
phiSpatialObserved(validPhi) = ...
    (xSymBulk(validPhi)-xA(validPhi))./den(validPhi);

zeroCrossingFraction = NaN(nTime,1);
for k = 1:nTime
    if JAS(k) ~= JSE(k)
        z0 = JAS(k)/(JAS(k)-JSE(k));
        if z0 >= 0 && z0 <= 1
            zeroCrossingFraction(k) = z0;
        end
    end
end

%% ========================================================================
% Outputs
% ========================================================================

out = struct;

out.integrationMode = 'intervalWiseExactD_spatialSymplast';
out.hydraulicFluxConvention = 'backwardFD_endingNode';
out.spatialSymplastCells = nCells;
out.spatialSymplastPathLength_m = pathLength;
out.spatialSymplastCellLength_m = dz;
out.spatialSymplastCellCenter_m = zCenter_m;
out.spatialSymplastFacePosition_m = zFace_m;

out.apoplasticWaterAmount_mol_m2 = WA.';
out.symplasticWaterAmount_mol_m2 = WS.';
out.symplasticCellWaterAmount_mol_m2 = Wcell;

out.JxylemToApoplast_mol_m2_s = Jx.';
out.JapoplastToSymplast_mol_m2_s = JAS.';
out.JsymplastToEvaporation_mol_m2_s = JSE.';
out.JapoplastToEvaporation_mol_m2_s = JAE.';
out.faceWaterFlux_mol_m2_s = faceWaterFlux;
out.zeroCrossingFraction = zeroCrossingFraction.';

out.ratioApoplast = rA.';
out.ratioSymplastCell = rCell;
out.ratioSymplastBulk = rSymBulk.';
out.ratioSymplast = rSymBulk.'; % compatibility alias
out.ratioEvaporationSite = rEvap.';

out.deltaDApoplast = deltaDA.';
out.deltaDSymplastCell = deltaDCell;
out.deltaDSymplastBulk = deltaDSymBulk.';
out.deltaDSymplast = deltaDSymBulk.'; % compatibility alias
out.deltaDEvaporationSite = deltaDEvap.';

out.DInventoryApoplast_molD_m2 = NDA.';
out.DInventorySymplastCell_molD_m2 = NDCell;
out.DInventorySymplast_molD_m2 = NDsym.';
out.DInventoryTotal_molD_m2 = NDtotal.';

out.cumulativeXylemDInput_molD_m2 = cumX.';
out.cumulativeAtmosphericDInput_molD_m2 = cumAtmIn.';
out.cumulativeAtmosphericDOutput_molD_m2 = cumAtmOut.';
out.cumulativeDirectApoplastEvapDTransfer_molD_m2 = cumAE.';
out.cumulativeSymplastFaceDTransfer_molD_m2 = cumFace;
out.faceDFlux_molD_m2_s = faceDFlux;
out.directApoplastEvapDFlux_molD_m2_s = directAEDFlux.';

out.DInventoryClosureResidual_molD_m2 = DInventoryClosureResidual.';
out.apoplastDClosureResidual_molD_m2 = apoplastDClosureResidual.';
out.symplastCellDClosureResidual_molD_m2 = cellDClosureResidual;
out.symplastDClosureResidual_molD_m2 = symplastDClosureResidual.';
out.evaporationSiteDClosureResidual_molD_m2_s = siteDClosure.';

out.apoplastWaterClosureResidual_mol_m2_s = apoplastWaterClosure.';
out.symplastWaterClosureResidual_mol_m2_s = symplastWaterClosure.';
out.wholeLeafWaterClosureResidual_mol_m2_s = wholeLeafWaterClosure.';
out.pathwayFlowResidual_mol_m2_s = pathwayFlowResidual.';
out.hydraulicPathwayResidual_mol_m2_s = hydraulicPathwayResidual.';
out.intervalApoplastForcingResidual_mol_m2_s = ...
    intervalApoplastForcingResidual.';
out.intervalSymplastForcingResidual_mol_m2_s = ...
    intervalSymplastForcingResidual.';
out.atmosphericWaterClosureResidual_mol_m2_s = atmosphericWaterClosure.';
out.xylemBoundaryFlowResidual_mol_m2_s = boundaryResidual.';

out.alphaEquilibrium = alphaEq.';
out.alphaKinetic = alphaK.';
out.ratioVaporAmbient = rV.';
out.ratioVaporIncomingEffective = rVin.';
out.atmosphericGrossInflow_mol_m2_s = Jin.';
out.atmosphericGrossOutflow_mol_m2_s = Jout.';

out.HDOdiffusivity_m2_s = DHDO.';
out.liquidWaterMolarConcentration_mol_m3 = Cwater.';
out.symplastDiffusiveFraction = fSymDiff.';
out.symplastWholePathDiffusionConductance_mol_m2_s = ...
    KdiffSymWholePath.';
out.pecletNumberSpatialBoundary = PeSymBoundary.';
out.pecletTransferFactorClassical = phiClassical.';
out.pecletTransferFactorSpatialObserved = phiSpatialObserved.';
KdiffWholeReference = Cwater.*DHDO.*aEffBase./pathLength;
PeWholeReference = E./KdiffWholeReference;
phiWholeReference = ones(nTime,1);
maskWhole = abs(PeWholeReference) > 1e-10;
phiWholeReference(maskWhole) = ...
    (1-exp(-PeWholeReference(maskWhole)))./PeWholeReference(maskWhole);
out.pecletNumberWholePathReference = PeWholeReference.';
out.pecletTransferFactorSteadyReference = phiWholeReference.';
out.leafAnatomyTransport = Aleaf;

out.odeRelTol = odeRelTol;
out.odeAbsTol = odeAbsTol;
out.odeMaxStepFraction = odeMaxStepFraction;
out.spatialSymplastModel = ...
    'finite-volume signed Scharfetter-Gummel D inventory pathway';
out.apoplastSymplastExchangeModel = ...
    'spatial advection-diffusion only; scalar GAS removed';

%% ========================================================================
% Nested RHS
% ========================================================================

    function dydt = spatialLeafRHS(tt,y,interval)

        f = (tt-interval.t0)/(interval.t1-interval.t0);
        f = min(max(f,0),1);

        WA_t = interval.WA0 + f*(interval.WA1-interval.WA0);
        WS_t = interval.WS0 + f*(interval.WS1-interval.WS0);
        Wcell_t = WS_t/interval.nCells;

        if WA_t <= 0 || Wcell_t <= 0
            error('Spatial leaf ODE encountered nonpositive water storage.');
        end

        NDA_t = y(interval.idxA);
        NDCell_t = y(interval.idxCell);

        xA_t = NDA_t/(2*WA_t);
        xCell_t = NDCell_t./(2*Wcell_t);

        if ~isfinite(xA_t) || xA_t <= 0 || xA_t >= 1 ...
                || any(~isfinite(xCell_t) | xCell_t <= 0 | xCell_t >= 1)
            error('Spatial leaf ODE produced an invalid D atom fraction.');
        end

        [~,xE_t] = solveSpatialEvaporationSiteExact( ...
            dAtomFractionToNormalizedRatio(xA_t,interval.Rstd), ...
            xCell_t(end), ...
            interval.rVin, ...
            interval.JAE, ...
            interval.Jin, ...
            interval.Jout, ...
            interval.KdiffApo, ...
            interval.Gright, ...
            interval.BminusRight, ...
            interval.BplusRight, ...
            interval.alphaEq, ...
            interval.alphaK, ...
            interval.Rstd);

        faceFlux_xD = spatialFaceFluxesSignedPrecomputed( ...
            xCell_t, ...
            xA_t, ...
            xE_t, ...
            interval.Jface, ...
            interval.Ginternal, ...
            interval.Gright, ...
            interval.BminusInternal, ...
            interval.BplusInternal, ...
            interval.BminusRight, ...
            interval.BplusRight);

        faceDFlux_t = 2.*faceFlux_xD;

        xX_t = ratioNormalizedToDAtomFraction(interval.rX,interval.Rstd);
        xVin_t = ratioNormalizedToDAtomFraction(interval.rVin,interval.Rstd);
        rE_t = dAtomFractionToNormalizedRatio(xE_t,interval.Rstd);
        rVout_t = rE_t/(interval.alphaEq*interval.alphaK);
        xVout_t = ratioNormalizedToDAtomFraction(rVout_t,interval.Rstd);

        FDx_t = 2*interval.Jx*xX_t;
        FDatmIn_t = 2*interval.Jin*xVin_t;
        FDatmOut_t = 2*interval.Jout*xVout_t;

        JAin_t = max(interval.JAE,0);
        JAout_t = max(-interval.JAE,0);
        FDAE_t = 2*( ...
            (JAin_t+interval.KdiffApo)*xA_t ...
            - (JAout_t+interval.KdiffApo)*xE_t);

        dNDA_t = FDx_t - faceDFlux_t(1) - FDAE_t;
        dNDCell_t = ...
            faceDFlux_t(1:interval.nCells) ...
            - faceDFlux_t(2:interval.nCells+1);

        dydt = zeros(size(y));
        dydt(interval.idxA) = dNDA_t;
        dydt(interval.idxCell) = dNDCell_t;
        dydt(interval.idxCumX) = FDx_t;
        dydt(interval.idxCumAtmIn) = FDatmIn_t;
        dydt(interval.idxCumAtmOut) = FDatmOut_t;
        dydt(interval.idxCumAE) = FDAE_t;
        dydt(interval.idxCumFace) = faceDFlux_t;

        % ode15s independent variable is days; physical fluxes are per s.
        dydt = interval.secondsPerDay.*dydt;
    end

end

%% ========================================================================
% Exact massless-site solution coupled to the spatial right boundary
% ========================================================================

function [rEvap,xE] = solveSpatialEvaporationSiteExact( ...
    rApo,xLast,rVin,JAE,Jin,Jout,KdiffApo, ...
    Gright,BminusRight,BplusRight,alphaEq,alphaK,Rstd)

xA = ratioNormalizedToDAtomFraction(rApo,Rstd);
xVin = ratioNormalizedToDAtomFraction(rVin,Rstd);
a = alphaEq*alphaK;

JAin = max(JAE,0);
JAout = max(-JAE,0);

% Site balance in xD-flux units:
%
%   (JAin+KA)*xA
% + Gright*B(-p)*xLast
% + Jin*xVin
%
% =
%
%   (JAout+KA+Gright*B(p))*xE
% + Jout*xVout.
%
% xE = q/(1+q), xVout = q/(a+q), q = absolute R_E.

B = ...
    (JAin+KdiffApo)*xA ...
    + Gright*BminusRight*xLast ...
    + Jin*xVin;

C = ...
    JAout ...
    + KdiffApo ...
    + Gright*BplusRight;

A2 = B-C-Jout;
A1 = B*(1+a)-C*a-Jout;
A0 = B*a;

scale = max([abs(A2),abs(A1),abs(A0),1]);

if abs(A2) < 1e-14*scale
    if abs(A1) < 1e-14*scale
        error('Degenerate spatial evaporation-site isotope balance.');
    end
    qCandidates = -A0/A1;
else
    disc = A1^2-4*A2*A0;
    if disc < -1e-14*scale^2
        error('Negative spatial evaporation-site quadratic discriminant.');
    end
    disc = max(disc,0);
    qCandidates = [ ...
        (-A1+sqrt(disc))/(2*A2), ...
        (-A1-sqrt(disc))/(2*A2)];
end

valid = isfinite(qCandidates) & isreal(qCandidates) & qCandidates > 0;
qCandidates = real(qCandidates(valid));
if isempty(qCandidates)
    error('No positive physical spatial evaporation-site isotope solution.');
end

residual = NaN(size(qCandidates));
for k = 1:numel(qCandidates)
    q = qCandidates(k);
    xEk = q/(1+q);
    xVout = q/(a+q);
    residual(k) = ...
        B ...
        - C*xEk ...
        - Jout*xVout;
end

[~,iBest] = min(abs(residual));
qEvap = qCandidates(iBest);
xE = qEvap/(1+qEvap);
rEvap = qEvap/Rstd;

end

%% ========================================================================
% Signed spatial face isotope fluxes
% ========================================================================

function faceFlux_xD = spatialFaceFluxesSignedPrecomputed( ...
    xCell,xA,xE,Jface,Ginternal,Gright, ...
    BminusInternal,BplusInternal,BminusRight,BplusRight)

nCells = numel(xCell);
faceFlux_xD = zeros(nCells+1,1);

% Left/apoplast boundary: Danckwerts total-flux inflow for JAS>=0;
% donor-aware outflow for JAS<0.
if Jface(1) >= 0
    faceFlux_xD(1) = Jface(1)*xA;
else
    faceFlux_xD(1) = Jface(1)*xCell(1);
end

for i = 1:(nCells-1)
    faceFlux_xD(i+1) = ...
        Ginternal*( ...
            BminusInternal(i)*xCell(i) ...
            - BplusInternal(i)*xCell(i+1));
end

faceFlux_xD(end) = ...
    Gright*( ...
        BminusRight*xCell(end) ...
        - BplusRight*xE);

end

function Jface = buildFaceWaterFlux(JAS,JSE,nCells)
fraction = (0:nCells).'./nCells;
Jface = JAS-fraction.*(JAS-JSE);
end

%% ========================================================================
% Ratio / forcing helpers
% ========================================================================

function x = ratioNormalizedToDAtomFraction(r,Rstd)
Rabs = r.*Rstd;
x = Rabs./(1+Rabs);
end

function r = dAtomFractionToNormalizedRatio(x,Rstd)
Rabs = x./(1-x);
r = Rabs./Rstd;
end

function x = expandTimeLocal(x,nTime)
if isscalar(x)
    x = repmat(x,nTime,1);
else
    x = x(:);
    if numel(x) ~= nTime
        error('Time-varying forcing does not match model time grid.');
    end
end
end

function x = expandScalarOrCellVector(x,nCells,name)
if isscalar(x)
    x = repmat(x,nCells,1);
else
    x = x(:);
    if numel(x) ~= nCells
        error('%s must be scalar or contain nCells values.',name);
    end
end
end

function value = getCfgScalar(S,name,defaultValue)
if isfield(S,name)
    value = S.(name);
else
    value = defaultValue;
end
if ~isscalar(value) || ~isfinite(value)
    error('Configuration field %s must be a finite scalar.',name);
end
end

function D = hdoLiquidDiffusivityLocal(TC)
TK = TC+273.15;
D = 116e-9.*exp(-626./(TK-139));
end

function es = saturationVaporPressureKPaLocal(TC)
es = 0.61078.*exp((17.2694.*TC)./(TC+237.29));
end

function B = bernoulliFunction(z)
B = zeros(size(z));
small = abs(z) < 1e-7;
zs = z(small);
B(small) = 1-zs./2+zs.^2./12-zs.^4./720;
largePos = z > 50;
B(largePos) = z(largePos).*exp(-z(largePos));
largeNeg = z < -50;
B(largeNeg) = -z(largeNeg);
mid = ~(small|largePos|largeNeg);
B(mid) = z(mid)./expm1(z(mid));
end
