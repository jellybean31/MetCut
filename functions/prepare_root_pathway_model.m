function P = prepare_root_pathway_model( ...
    rootTraits, ...
    soilHydro, ...
    plantState, ...
    psiMinMPa, ...
    nLayer)
%PREPARE_ROOT_PATHWAY_MODEL
%
% Precompute static root geometry, radial root conductance,
% matric-flux-potential tables, and the monotonic H(psi) lookup used
% to solve root-surface potentials efficiently.

%% Root traits

RLD = expandLayer( ...
    rootTraits.rootLengthDensity_m_m3, ...
    nLayer, ...
    'rootLengthDensity_m_m3');

rRoot = expandLayer( ...
    rootTraits.rootRadius_m, ...
    nLayer, ...
    'rootRadius_m');

Lp = expandLayer( ...
    rootTraits.Lp_m_s_MPa, ...
    nLayer, ...
    'Lp_m_s_MPa');

dz = expandLayer( ...
    soilHydro.layerThickness_m, ...
    nLayer, ...
    'layerThickness_m');

%% Optional sensitivity scaling

if isfield(rootTraits,'RLDScale')
    RLD = RLD .* rootTraits.RLDScale;
end

if isfield(rootTraits,'rootRadiusScale')
    rRoot = rRoot .* rootTraits.rootRadiusScale;
end

if isfield(rootTraits,'LpScale')
    Lp = Lp .* rootTraits.LpScale;
end

%% LAI

if ~isfield(plantState,'leafAreaIndex')
    error('plantState.leafAreaIndex is required.');
end

%% Hydraulic redistribution

if isfield(rootTraits,'allowHydraulicRedistribution')

    if ~isscalar(rootTraits.allowHydraulicRedistribution)

        error( ...
            'rootTraits.allowHydraulicRedistribution must be scalar.');

    end

    allowHydraulicRedistribution = ...
        logical(rootTraits.allowHydraulicRedistribution);

else

    % Historical uptake-only behavior.

    allowHydraulicRedistribution = ...
        false;

end

%% Optional hydraulic-redistribution capacity scaling
% Default = 1 preserves established production behavior.
if isfield(rootTraits,'hydraulicRedistributionScale')
    hydraulicRedistributionScale = rootTraits.hydraulicRedistributionScale;
    if ~isscalar(hydraulicRedistributionScale) ...
            || ~isfinite(hydraulicRedistributionScale) ...
            || hydraulicRedistributionScale < 0
        error('rootTraits.hydraulicRedistributionScale must be a finite nonnegative scalar.');
    end
else
    hydraulicRedistributionScale = 1;
end

%% Optional diagnostic lower bracket for nonlinear root-demand solve
% Default = -5 MPa preserves established production behavior.
if isfield(rootTraits,'rootDemandPsiMinMPa')
    rootDemandPsiMinMPa = rootTraits.rootDemandPsiMinMPa;
    if ~isscalar(rootDemandPsiMinMPa) ...
            || ~isfinite(rootDemandPsiMinMPa) ...
            || rootDemandPsiMinMPa >= 0
        error('rootTraits.rootDemandPsiMinMPa must be a finite negative scalar.');
    end
else
    rootDemandPsiMinMPa = -5;
end

LAI = plantState.leafAreaIndex;

%% Constants

rhoMolWater = 55500;    % mol m^-3

rhoWater = 1000;        % kg m^-3
g = 9.80665;            % m s^-2

headPerMPa = ...
    1e6 / (rhoWater*g); % m MPa^-1

%% Root geometry

rootLengthPerGroundArea = ...
    RLD .* dz;

rootLengthPerLeafArea = ...
    rootLengthPerGroundArea ./ LAI;

rootSurfaceAreaPerLeafArea = ...
    2 .* pi .* rRoot .* ...
    rootLengthPerLeafArea;

geometryRootFraction = ...
    rootLengthPerGroundArea ...
    ./ sum(rootLengthPerGroundArea);

%% Root tissue conductance

KrootTissue = ...
    Lp .* ...
    rootSurfaceAreaPerLeafArea .* ...
    rhoMolWater;

%% Cylindrical soil domain

rOuter = ...
    sqrt( ...
        1 ./ (pi .* RLD));

rOuter = ...
    max( ...
        rOuter, ...
        1.01 .* rRoot);

logGeometry = ...
    log(rOuter ./ rRoot);

%% Matric flux potential

F = soil_matric_flux_potential( ...
    soilHydro, ...
    psiMinMPa, ...
    nLayer);

%% Radial soil-flow geometry coefficient
%
% Soil molar flux:
%
%   Jsoil =
%       Csoil * [Phi(psiSoil)-Phi(psiRootSurface)]

Csoil = ...
    2 .* pi .* ...
    rootLengthPerLeafArea .* ...
    rhoMolWater .* ...
    headPerMPa ...
    ./ logGeometry;

%% Monotonic H-function
%
% Equality of soil and root tissue flux gives:
%
% Kroot*(psiRS-psiRoot)
%   =
% Csoil*(PhiSoil-PhiRS)
%
% therefore:
%
% Kroot*psiRS + Csoil*PhiRS
%   =
% Kroot*psiRoot + Csoil*PhiSoil
%
% Define the left side as H(psiRS).

psiGrid = F.psiGrid;

Hgrid = NaN( ...
    nLayer, ...
    numel(psiGrid));

for i = 1:nLayer

    Hgrid(i,:) = ...
        KrootTissue(i) .* psiGrid ...
        + ...
        Csoil(i) .* F.Phi(i,:);

end

%% Output

P.nLayer = nLayer;

P.RLD = RLD;
P.rootRadius_m = rRoot;
P.Lp_m_s_MPa = Lp;

P.rootLengthPerGroundArea = ...
    rootLengthPerGroundArea;

P.rootLengthPerLeafArea = ...
    rootLengthPerLeafArea;

P.rootSurfaceAreaPerLeafArea = ...
    rootSurfaceAreaPerLeafArea;

P.geometryRootFraction = ...
    geometryRootFraction;

P.rOuter_m = rOuter;
P.logGeometry = logGeometry;

P.KrootTissue = KrootTissue;

P.Csoil = Csoil;

P.mfp = F;

P.Hgrid = Hgrid;

P.allowHydraulicRedistribution = ...
    allowHydraulicRedistribution;

P.hydraulicRedistributionScale = ...
    hydraulicRedistributionScale;

P.rootDemandPsiMinMPa = ...
    rootDemandPsiMinMPa;

end


%% ========================================================================
function x = expandLayer(x,nLayer,name)

if isscalar(x)
    x = repmat(x,nLayer,1);
else
    x = x(:);
end

if numel(x) ~= nLayer

    error( ...
        '%s must be scalar or contain one value per soil layer.', ...
        name);

end

end