function L = soil_root_layer_flux( ...
    psiSoil, ...
    psiRoot, ...
    P)
%SOIL_ROOT_LAYER_FLUX
%
% Compute nonlinear radial soil -> root flux for every soil layer at a
% specified root-collar potential.
%
% Reverse flow is disabled.
%
% INPUT
%
%   psiSoil  = nLayer x 1 bulk soil potentials, MPa
%   psiRoot  = root-collar potential, MPa
%   P        = output from prepare_root_pathway_model
%
% OUTPUT
%
%   flux
%   psiRootSurface
%   KsoilEquivalent
%   KrootTissue
%   Keffective
%   soilLimitationFraction
%   rootLimitationFraction
%   relativeSoilConductivity
%   KsoilUnsat_m_s
%   KsoilPathBulk_m_s
%   soilFluxResidual
%   active

allowHR = ...
    isfield(P,'allowHydraulicRedistribution') ...
    && logical(P.allowHydraulicRedistribution);

if isfield(P,'hydraulicRedistributionScale')
    redistributionScale = P.hydraulicRedistributionScale;
else
    redistributionScale = 1;
end

psiSoil = psiSoil(:);

nLayer = P.nLayer;

isUptake = ...
    false(nLayer,1);


isRelease = ...
    false(nLayer,1);


flowDirection = ...
    zeros(nLayer,1);

%% Preallocate

flux = zeros(nLayer,1);

psiRS = NaN(nLayer,1);

KsoilEq = NaN(nLayer,1);

Keff = zeros(nLayer,1);

soilLim = NaN(nLayer,1);
rootLim = NaN(nLayer,1);

KrelBulk = NaN(nLayer,1);
KsoilUnsatBulk = NaN(nLayer,1);
KsoilPathBulk = NaN(nLayer,1);

fluxResidual = NaN(nLayer,1);

active = false(nLayer,1);

psiGrid = P.mfp.psiGrid;

tolPsi = 1e-12;

%% ------------------------------------------------------------------------
% Layers
% -------------------------------------------------------------------------

for i = 1:nLayer

    psiS = psiSoil(i);

    if ~isfinite(psiS) || ~isfinite(psiRoot)
        continue
    end

    %% Bulk-soil diagnostics

    KrelBulk(i) = interp1( ...
        psiGrid, ...
        P.mfp.Krel(i,:), ...
        psiS, ...
        'linear');

    KsoilUnsatBulk(i) = interp1( ...
        psiGrid, ...
        P.mfp.KsoilUnsat_m_s(i,:), ...
        psiS, ...
        'linear');

    KsoilPathBulk(i) = interp1( ...
        psiGrid, ...
        P.mfp.KsoilPath_m_s(i,:), ...
        psiS, ...
        'linear');

    %% ================================================================
    % Flow-direction control
    % ================================================================

    if ...
            ~allowHR ...
            && psiS <= psiRoot

        % Historical uptake-only behavior.

        psiRS(i) = ...
            psiS;


        KsoilEq(i) = ...
            P.Csoil(i) ...
            .* KsoilPathBulk(i);


        continue

    end

    %% Matric flux potential at bulk soil

    PhiS = interp1( ...
        psiGrid, ...
        P.mfp.Phi(i,:), ...
        psiS, ...
        'linear');

    %% Solve root-surface potential by inverting H
    %
    % H(psiRS) =
    %   Kroot*psiRoot + Csoil*Phi(psiSoil)

    rhs = ...
        P.KrootTissue(i) .* psiRoot ...
        + ...
        P.Csoil(i) .* PhiS;

    H = P.Hgrid(i,:);

    % Numerical protection

    rhs = min(max(rhs,H(1)),H(end));

    psiRS(i) = interp1( ...
        H, ...
        psiGrid, ...
        rhs, ...
        'linear');

    % Root surface must remain between root xylem and bulk soil.

    psiLo = ...
        min(psiRoot,psiS);


    psiHi = ...
        max(psiRoot,psiS);


    psiRS(i) = ...
        min( ...
        max(psiRS(i),psiLo), ...
        psiHi);

    %% Root-tissue flux

    Jroot = ...
        P.KrootTissue(i) ...
        .* (psiRS(i)-psiRoot);


    if ~allowHR

        Jroot = ...
            max(Jroot,0);

    end

    %% Soil flux independently calculated for mass-balance diagnostic

    PhiRS = interp1( ...
        psiGrid, ...
        P.mfp.Phi(i,:), ...
        psiRS(i), ...
        'linear');

    Jsoil = ...
        P.Csoil(i) ...
        .* (PhiS-PhiRS);


    if ~allowHR

        Jsoil = ...
            max(Jsoil,0);

    elseif Jroot < 0

        % Scale only reverse root -> soil transfer capacity. Forward uptake
        % remains unchanged. Apply the same factor to both pathway sides so
        % the layer soil/root flux residual remains a closure diagnostic.
        Jroot = redistributionScale .* Jroot;
        Jsoil = redistributionScale .* Jsoil;

    end

    %% Use root-side flux

    flux(i) = Jroot;

    flowTol = ...
        1e-12;


    if flux(i) > flowTol

        isUptake(i) = ...
            true;


        flowDirection(i) = ...
            1;


    elseif flux(i) < -flowTol

        isRelease(i) = ...
            true;


        flowDirection(i) = ...
            -1;

    end


    active(i) = ...
        abs(flux(i)) > flowTol;

    fluxResidual(i) = ...
        Jsoil-Jroot;

    %% Potential drops

    soilDrop = ...
        psiS-psiRS(i);

    rootDrop = ...
        psiRS(i)-psiRoot;

    totalDrop = ...
        psiS-psiRoot;

    %% Secant-equivalent conductances

    if abs(soilDrop) > tolPsi

        KsoilEq(i) = ...
            flux(i) ./ soilDrop;

    else

        KsoilEq(i) = ...
            P.Csoil(i) ...
            .* KsoilPathBulk(i);

    end


    if abs(totalDrop) > tolPsi

        Keff(i) = ...
            flux(i) ./ totalDrop;


        soilLim(i) = ...
            soilDrop ./ totalDrop;


        rootLim(i) = ...
            rootDrop ./ totalDrop;

    else

        Keff(i) = ...
            0;

    end

end

%% Output

L.flux = flux;

L.psiRootSurface = psiRS;

L.KsoilEquivalent = KsoilEq;

L.KrootTissue = P.KrootTissue;

L.Keffective = Keff;

L.soilLimitationFraction = soilLim;
L.rootLimitationFraction = rootLim;

L.relativeSoilConductivity = KrelBulk;

L.KsoilUnsat_m_s = KsoilUnsatBulk;

L.KsoilPathBulk_m_s = ...
    KsoilPathBulk;

L.soilFluxResidual = ...
    fluxResidual;

L.isUptake = ...
    isUptake;


L.isRelease = ...
    isRelease;


L.flowDirection = ...
    flowDirection;


L.allowHydraulicRedistribution = ...
    allowHR;

L.totalUptakeFlux = ...
    sum( ...
        max(flux,0), ...
        'omitnan');


L.totalReleaseFlux = ...
    sum( ...
        max(-flux,0), ...
        'omitnan');


L.netRootFlux = ...
    sum( ...
        flux, ...
        'omitnan');

L.active = active;

end