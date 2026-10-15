function F = soil_matric_flux_potential( ...
    soilHydro, ...
    psiMinMPa, ...
    nLayer)
%SOIL_MATRIC_FLUX_POTENTIAL
%
% Precompute soil hydraulic conductivity and matric flux potential:
%
%   Phi(psi) = integral K(psi) dpsi
%
% Phi is evaluated in pressure-potential space (MPa). Conversion from
% MPa to hydraulic head is applied later in the radial-flow geometry.
%
% OUTPUT
%
%   F.psiGrid
%   F.Krel
%   F.KsoilUnsat_m_s
%   F.KsoilPath_m_s
%   F.Phi
%
% Phi units:
%   m s^-1 MPa

%% Grid size

if isfield(soilHydro,'mfpGridSize')
    nGrid = soilHydro.mfpGridSize;
else
    nGrid = 4000;
end

nGrid = max(round(nGrid),500);

%% Layer parameters

Ksat = expandLayer( ...
    soilHydro.Ksat_m_s, ...
    nLayer, ...
    'Ksat_m_s');

alpha = expandLayer( ...
    soilHydro.alphaMPa, ...
    nLayer, ...
    'alphaMPa');

nVG = expandLayer( ...
    soilHydro.n, ...
    nLayer, ...
    'n');

lVG = expandLayer( ...
    soilHydro.l, ...
    nLayer, ...
    'l');

if isfield(soilHydro,'rhizosphereFactor')

    rhizo = expandLayer( ...
        soilHydro.rhizosphereFactor, ...
        nLayer, ...
        'rhizosphereFactor');

else

    rhizo = ones(nLayer,1);

end

if isfield(soilHydro,'KsoilScale')

    scale = expandLayer( ...
        soilHydro.KsoilScale, ...
        nLayer, ...
        'KsoilScale');

else

    scale = ones(nLayer,1);

end

%% Validation

if any(Ksat <= 0)
    error('Ksat_m_s must be positive.');
end

if any(alpha <= 0)
    error('alphaMPa must be positive.');
end

if any(nVG <= 1)
    error('van Genuchten n must be > 1.');
end

%% Potential grid
%
% Log spacing in suction provides much better resolution near saturation
% than a uniform psi grid.

psiTiny = 1e-8;

psiMinMPa = min(psiMinMPa,-10*psiTiny);

suction = logspace( ...
    log10(abs(psiMinMPa)), ...
    log10(psiTiny), ...
    nGrid-1);

psiGrid = [ ...
    -suction, ...
    0];

%% Preallocate

Krel = NaN(nLayer,nGrid);

KsoilUnsat = NaN(nLayer,nGrid);

KsoilPath = NaN(nLayer,nGrid);

Phi = NaN(nLayer,nGrid);

%% Build curves

for i = 1:nLayer

    Krel(i,:) = ...
        soil_relative_conductivity( ...
            psiGrid, ...
            alpha(i), ...
            nVG(i), ...
            lVG(i));

    % Bulk soil conductivity

    KsoilUnsat(i,:) = ...
        Ksat(i) .* Krel(i,:);

    % Conductivity actually used by soil/root pathway

    KsoilPath(i,:) = ...
        KsoilUnsat(i,:) ...
        .* rhizo(i) ...
        .* scale(i);

    % Matric flux potential

    Phi(i,:) = ...
        cumtrapz( ...
            psiGrid, ...
            KsoilPath(i,:));

end

%% Output

F.psiGrid = psiGrid;

F.Krel = Krel;

F.KsoilUnsat_m_s = ...
    KsoilUnsat;

F.KsoilPath_m_s = ...
    KsoilPath;

F.Phi = Phi;

F.Ksat_m_s = Ksat;
F.alphaMPa = alpha;
F.n = nVG;
F.l = lVG;

F.rhizosphereFactor = rhizo;
F.KsoilScale = scale;

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