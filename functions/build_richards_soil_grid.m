function [G,soilGrid,rootGrid,psiInitialGrid_MPa] = ...
    build_richards_soil_grid( ...
        soilHydro, ...
        rootTraits, ...
        psiInitialZone_MPa, ...
        edges_m)
%BUILD_RICHARDS_SOIL_GRID
%
% Construct a vertically resolved 1-D soil grid while preserving the
% properties of the existing three soil/root zones.
%
% Each hydraulic cell inherits the constitutive and root properties of
% its parent zone.
%
% rootLengthDensity is not divided among cells.
%
% prepare_root_pathway_model later calculates:
%
%   root length per ground area = RLD * cell thickness
%
% Therefore subdivision automatically preserves integrated root length:
%
%   sum_i RLD_i * dz_i
%
% provided the original RLD is interpreted as a volumetric density.
%
%
% INPUT
%
% soilHydro
%   Existing three-zone soil hydraulic parameter structure.
%
% rootTraits
%   Existing three-zone root parameter structure.
%
% psiInitialZone_MPa
%   Three initial soil water potentials.
%
% edges_m
%   Optional vector of cell edges from soil surface downward.
%
% If omitted, a graded grid is generated automatically with high
% resolution near the surface.
%
%
% OUTPUT
%
% G
%   Grid geometry and parent-zone mapping.
%
% soilGrid
%   soilHydro expanded to one value per hydraulic cell.
%
% rootGrid
%   rootTraits expanded to one value per hydraulic cell.
%
% psiInitialGrid_MPa
%   Initial water potential for each hydraulic cell.


%% ========================================================================
% Original geometry
% ========================================================================

zOriginal = ...
    soilHydro.layerThickness_m(:);

if numel(zOriginal) ~= 3

    error( ...
        'Current grid builder expects three original soil zones.');

end

if any(~isfinite(zOriginal)) || any(zOriginal <= 0)

    error( ...
        'soilHydro.layerThickness_m must contain positive values.');

end


zoneBottom_m = ...
    cumsum(zOriginal);

zTotal_m = ...
    zoneBottom_m(end);


%% ========================================================================
% Initial potentials
% ========================================================================

psiInitialZone = ...
    psiInitialZone_MPa(:);

if numel(psiInitialZone) ~= 3

    error( ...
        'psiInitialZone_MPa must contain three values.');

end

if any(~isfinite(psiInitialZone))

    error( ...
        'psiInitialZone_MPa must contain finite values.');

end


%% ========================================================================
% Default graded grid
% ========================================================================

if nargin < 4 || isempty(edges_m)

    z1 = ...
        zoneBottom_m(1);

    z2 = ...
        zoneBottom_m(2);


    % High resolution is concentrated near the atmosphere/soil
    % interface, where rainfall wetting fronts and isotope gradients
    % should be strongest.
    %
    % Current 30-cm first-layer case gives approximately:
    %
    %   0-5 cm      : 0.5-cm cells
    %   5-15 cm     : 1-cm cells
    %   15-30 cm    : 2.5-cm cells
    %   30-60 cm    : 5-cm cells
    %   60-90 cm    : 5-cm cells

    topFineEnd = ...
        min(0.05,z1);

    topMediumEnd = ...
        min(0.15,z1);


    edges1a = ...
        makeSegmentEdges( ...
            0, ...
            topFineEnd, ...
            0.005);


    edges1b = ...
        makeSegmentEdges( ...
            topFineEnd, ...
            topMediumEnd, ...
            0.010);


    edges1c = ...
        makeSegmentEdges( ...
            topMediumEnd, ...
            z1, ...
            0.025);


    edges2 = ...
        makeSegmentEdges( ...
            z1, ...
            z2, ...
            0.050);


    edges3 = ...
        makeSegmentEdges( ...
            z2, ...
            zTotal_m, ...
            0.050);


    edges_m = ...
        unique([ ...
            edges1a, ...
            edges1b, ...
            edges1c, ...
            edges2, ...
            edges3, ...
            0, ...
            zoneBottom_m.'], ...
            'sorted');

else

    edges_m = ...
        edges_m(:).';

end


%% ========================================================================
% Clean grid edges
%
% Floating-point colon operations can produce nearly duplicate values at
% segment boundaries (e.g. 0.05 versus 0.05000000000000001).
%
% Exact UNIQUE is insufficient because those values are not bitwise equal.
% Merge edges separated by less than edgeMergeTolerance.
% ========================================================================

edgeMergeTolerance = ...
    1e-10;


%% ------------------------------------------------------------------------
% Assemble requested edges plus mandatory domain/horizon boundaries
% -------------------------------------------------------------------------

edges_m = ...
    sort([ ...
        edges_m(:); ...
        0; ...
        zoneBottom_m(:)]);


%% ------------------------------------------------------------------------
% Basic validation
% -------------------------------------------------------------------------

if any(~isfinite(edges_m))

    error( ...
        'Grid edges must contain only finite values.');

end


if ...
        any(edges_m < -edgeMergeTolerance) ...
        || any(edges_m > zTotal_m + edgeMergeTolerance)

    error( ...
        'Grid edges must lie within the soil domain.');

end


%% ------------------------------------------------------------------------
% Snap values already numerically equal to domain boundaries
% -------------------------------------------------------------------------

edges_m( ...
    abs(edges_m) <= edgeMergeTolerance) = ...
    0;


edges_m( ...
    abs(edges_m - zTotal_m) <= edgeMergeTolerance) = ...
    zTotal_m;


%% ------------------------------------------------------------------------
% Merge near-duplicate edges
% -------------------------------------------------------------------------

keep = ...
    true(size(edges_m));


lastKept = ...
    1;


for i = 2:numel(edges_m)

    if ...
            edges_m(i) ...
            - edges_m(lastKept) ...
            <= edgeMergeTolerance

        keep(i) = ...
            false;

    else

        lastKept = ...
            i;

    end

end


edges_m = ...
    edges_m(keep);


%% ------------------------------------------------------------------------
% Force original horizon boundaries to their exact values
% -------------------------------------------------------------------------

for j = 1:numel(zoneBottom_m)

    [distance,iClosest] = ...
        min(abs( ...
            edges_m ...
            - zoneBottom_m(j)));


    if distance <= edgeMergeTolerance

        edges_m(iClosest) = ...
            zoneBottom_m(j);

    else

        edges_m = ...
            sort([ ...
                edges_m; ...
                zoneBottom_m(j)]);

    end

end


%% ------------------------------------------------------------------------
% Exact domain boundaries
% -------------------------------------------------------------------------

edges_m(1) = ...
    0;


edges_m(end) = ...
    zTotal_m;


%% ------------------------------------------------------------------------
% Final safety checks
% -------------------------------------------------------------------------

dzCheck = ...
    diff(edges_m);


if any(dzCheck <= edgeMergeTolerance)

    error( ...
        'Grid contains zero-width or near-zero-width cells.');

end


if ...
        abs(edges_m(1)) > edgeMergeTolerance ...
        || abs(edges_m(end) - zTotal_m) > edgeMergeTolerance

    error( ...
        'Grid does not span the complete soil domain.');

end

%% ========================================================================
% Cell geometry
% ========================================================================

zTop_m = ...
    edges_m(1:end-1);

zBottom_m = ...
    edges_m(2:end);

dz_m = ...
    zBottom_m ...
    - zTop_m;

zCenter_m = ...
    0.5 ...
    .* ( ...
        zTop_m ...
        + zBottom_m);


nCell = ...
    numel(dz_m);


%% ========================================================================
% Parent-zone assignment
% ========================================================================

zoneIndex = ...
    NaN(nCell,1);


for i = 1:nCell

    zc = ...
        zCenter_m(i);


    if zc <= zoneBottom_m(1)

        zoneIndex(i) = ...
            1;

    elseif zc <= zoneBottom_m(2)

        zoneIndex(i) = ...
            2;

    else

        zoneIndex(i) = ...
            3;

    end

end


%% ========================================================================
% Expanded soil hydraulic structure
% ========================================================================

soilGrid = ...
    soilHydro;


soilGrid.layerThickness_m = ...
    dz_m(:);


soilFields = { ...
    'thetaR', ...
    'thetaS', ...
    'alphaMPa', ...
    'n', ...
    'm', ...
    'l', ...
    'Ksat_m_s', ...
    'rhizosphereFactor', ...
    'KsoilScale'};


for q = 1:numel(soilFields)

    name = ...
        soilFields{q};


    if isfield(soilHydro,name)

        soilGrid.(name) = ...
            expandZoneToGrid( ...
                soilHydro.(name), ...
                zoneIndex, ...
                name);

    end

end


%% ========================================================================
% Expanded root structure
% ========================================================================

rootGrid = ...
    rootTraits;


rootFields = { ...
    'rootLengthDensity_m_m3', ...
    'rootRadius_m', ...
    'Lp_m_s_MPa', ...
    'RLDScale', ...
    'rootRadiusScale', ...
    'LpScale'};


for q = 1:numel(rootFields)

    name = ...
        rootFields{q};


    if isfield(rootTraits,name)

        value = ...
            rootTraits.(name);


        if isscalar(value)

            % Scalar scaling factors can remain scalar.
            %
            % Geometric/physiological fields are expanded so that the
            % resulting structure explicitly describes every cell.

            if any(strcmp( ...
                    name, ...
                    { ...
                    'rootLengthDensity_m_m3', ...
                    'rootRadius_m', ...
                    'Lp_m_s_MPa'}))

                rootGrid.(name) = ...
                    repmat( ...
                        value, ...
                        nCell, ...
                        1);

            else

                rootGrid.(name) = ...
                    value;

            end

        else

            rootGrid.(name) = ...
                expandZoneToGrid( ...
                    value, ...
                    zoneIndex, ...
                    name);

        end

    end

end


%% ========================================================================
% Initial cell potentials
% ========================================================================

psiInitialGrid_MPa = ...
    psiInitialZone(zoneIndex);


%% ========================================================================
% Integrated-root-geometry regression diagnostics
% ========================================================================

RLDoriginal = ...
    expandThree( ...
        rootTraits.rootLengthDensity_m_m3, ...
        'rootLengthDensity_m_m3');


RLDgrid = ...
    rootGrid.rootLengthDensity_m_m3(:);


% Include RLDScale because prepare_root_pathway_model applies this before
% calculating root length.

if isfield(rootTraits,'RLDScale')

    scaleOriginal = ...
        expandThree( ...
            rootTraits.RLDScale, ...
            'RLDScale');

else

    scaleOriginal = ...
        ones(3,1);

end


if isfield(rootGrid,'RLDScale')

    if isscalar(rootGrid.RLDScale)

        scaleGrid = ...
            repmat( ...
                rootGrid.RLDScale, ...
                nCell, ...
                1);

    else

        scaleGrid = ...
            rootGrid.RLDScale(:);

    end

else

    scaleGrid = ...
        ones(nCell,1);

end


rootLengthOriginal_m_m2 = ...
    RLDoriginal ...
    .* scaleOriginal ...
    .* zOriginal;


rootLengthGridByCell_m_m2 = ...
    RLDgrid ...
    .* scaleGrid ...
    .* dz_m(:);


rootLengthGridByZone_m_m2 = ...
    zeros(3,1);


for j = 1:3

    rootLengthGridByZone_m_m2(j) = ...
        sum( ...
            rootLengthGridByCell_m_m2( ...
                zoneIndex == j));

end


rootLengthRelativeError = ...
    ( ...
        rootLengthGridByZone_m_m2 ...
        - rootLengthOriginal_m_m2) ...
    ./ max( ...
        abs(rootLengthOriginal_m_m2), ...
        eps);


%% ========================================================================
% Output grid structure
% ========================================================================

G.nCell = ...
    nCell;

G.edges_m = ...
    edges_m(:);

G.zTop_m = ...
    zTop_m(:);

G.zBottom_m = ...
    zBottom_m(:);

G.zCenter_m = ...
    zCenter_m(:);

G.dz_m = ...
    dz_m(:);

G.totalDepth_m = ...
    zTotal_m;

G.originalLayerThickness_m = ...
    zOriginal;

G.originalZoneBottom_m = ...
    zoneBottom_m;

G.zoneIndex = ...
    zoneIndex;

G.rootLengthOriginal_m_m2 = ...
    rootLengthOriginal_m_m2;

G.rootLengthGridByZone_m_m2 = ...
    rootLengthGridByZone_m_m2;

G.rootLengthRelativeError = ...
    rootLengthRelativeError;

end


%% ========================================================================
function xGrid = expandZoneToGrid(x,zoneIndex,name)

x = ...
    x(:);


if isscalar(x)

    x = ...
        repmat(x,3,1);

end


if numel(x) ~= 3

    error( ...
        '%s must be scalar or contain one value per original zone.', ...
        name);

end


xGrid = ...
    x(zoneIndex);

end


%% ========================================================================
function x = expandThree(x,name)

x = ...
    x(:);


if isscalar(x)

    x = ...
        repmat(x,3,1);

end


if numel(x) ~= 3

    error( ...
        '%s must be scalar or contain three values.', ...
        name);

end

end


%% ========================================================================
function e = makeSegmentEdges(a,b,dz)

if b <= a

    e = ...
        a;

    return

end


e = ...
    a:dz:b;


if isempty(e) || abs(e(end)-b) > 1e-12

    e = ...
        [e b];

else

    e(end) = ...
        b;

end

end




