function S = leaf_storage_pressure_volume_state(psiLeaf,cfg)
%LEAF_STORAGE_PRESSURE_VOLUME_STATE
%
% Pressure-volume representation of leaf hydraulic water storage.
%
% INPUT
%   psiLeaf    leaf water potential [MPa]
%
%   cfg.saturatedWaterContent_mol_m2
%   cfg.osmoticPotentialFullTurgor_MPa
%   cfg.turgorLossPoint_MPa
%   cfg.relativeWaterContentTurgorLoss
%   cfg.apoplasticFraction
%
% Optional:
%   cfg.minimumRelativeWaterContent
%
%
% OUTPUT
%   S.amount_mol_m2
%   S.rawAmount_mol_m2
%   S.capacitanceEffective_mol_m2_MPa
%   S.boundActive
%
% Additional physiological diagnostics:
%   S.relativeWaterContent
%   S.symplasticRelativeWaterContent
%   S.apoplasticWaterAmount_mol_m2
%   S.symplasticWaterAmount_mol_m2
%   S.osmoticPotential_MPa
%   S.turgorPotential_MPa
%   S.reconstructedPsiLeaf_MPa
%   S.constitutiveResidual_MPa
%   S.belowTurgorLoss
%
% The formulation follows standard leaf pressure-volume
% relationships expressed in terms of total RWC.


%% ========================================================================
% Inputs
% ========================================================================

psi = ...
    psiLeaf;


Wsat = ...
    cfg.saturatedWaterContent_mol_m2;


psiO = ...
    cfg.osmoticPotentialFullTurgor_MPa;


psiTLP = ...
    cfg.turgorLossPoint_MPa;


rTLP = ...
    cfg.relativeWaterContentTurgorLoss;


af = ...
    cfg.apoplasticFraction;


%% ========================================================================
% Validation
% ========================================================================

if ...
        ~isscalar(Wsat) ...
        || ~isfinite(Wsat) ...
        || Wsat <= 0

    error( ...
        'saturatedWaterContent_mol_m2 must be positive and finite.');

end


if ...
        ~isscalar(psiO) ...
        || ~isfinite(psiO) ...
        || psiO >= 0

    error( ...
        'osmoticPotentialFullTurgor_MPa must be negative.');

end


if ...
        ~isscalar(psiTLP) ...
        || ~isfinite(psiTLP) ...
        || psiTLP >= psiO

    error([ ...
        'turgorLossPoint_MPa must be more negative than ' ...
        'osmoticPotentialFullTurgor_MPa.']);

end


if ...
        ~isscalar(rTLP) ...
        || ~isfinite(rTLP) ...
        || rTLP <= 0 ...
        || rTLP >= 1

    error( ...
        'relativeWaterContentTurgorLoss must lie in (0,1).');

end


if ...
        ~isscalar(af) ...
        || ~isfinite(af) ...
        || af < 0 ...
        || af >= rTLP

    error([ ...
        'apoplasticFraction must be nonnegative and less ' ...
        'than relativeWaterContentTurgorLoss.']);

end


if isfield(cfg,'minimumRelativeWaterContent')

    rMin = ...
        cfg.minimumRelativeWaterContent;

else

    % Prevent negative symplastic water.
    rMin = ...
        af ...
        + 1e-6;

end


if ...
        ~isscalar(rMin) ...
        || ~isfinite(rMin) ...
        || rMin <= af ...
        || rMin >= rTLP

    error([ ...
        'minimumRelativeWaterContent must be greater than ' ...
        'apoplasticFraction and less than RWC at turgor loss.']);

end


%% ========================================================================
% Reusable P-V constants
% ========================================================================

A = ...
    psiO ...
    .* psiTLP ...
    .* (1 - rTLP);


d0 = ...
    psiO ...
    - psiTLP .* rTLP;


d1 = ...
    psiTLP ...
    - psiO;


% Absolute slope of the linear turgor-RWC relation above TLP.

kTurgor = ...
    -psiO ...
    ./ (1 - rTLP);


%% ========================================================================
% Invert psiLeaf -> total RWC
% ========================================================================

rwcRaw = ...
    NaN(size(psi));


%% ------------------------------------------------------------------------
% Saturated / positive-pressure condition
% -------------------------------------------------------------------------

idxSat = ...
    psi >= 0;


rwcRaw(idxSat) = ...
    1;


%% ------------------------------------------------------------------------
% Below turgor loss
%
% psiLeaf = psiOsmotic because turgor = 0.
%
% This branch has an explicit linear-rational inverse.
% -------------------------------------------------------------------------

idxBelow = ...
    psi <= psiTLP;


idxBelow = ...
    idxBelow ...
    & ~idxSat;


rwcRaw(idxBelow) = ...
    ( ...
        A ./ psi(idxBelow) ...
        - d0 ...
    ) ...
    ./ d1;


%% ------------------------------------------------------------------------
% Above turgor loss
%
% psiLeaf =
%
%   A / (d0 + d1*RWC)
%   + kTurgor*(RWC - RWC_TLP)
%
% Rearrangement gives a quadratic in RWC.
% -------------------------------------------------------------------------

idxTurgid = ...
    psi < 0 ...
    & psi > psiTLP;


if any(idxTurgid,'all')

    psiNow = ...
        psi(idxTurgid);


    B = ...
        psiNow ...
        + kTurgor .* rTLP;


    q2 = ...
        -kTurgor .* d1;


    q1 = ...
        B .* d1 ...
        - kTurgor .* d0;


    q0 = ...
        B .* d0 ...
        - A;


    discriminant = ...
        q1.^2 ...
        - 4 .* q2 .* q0;


    if any(discriminant < -1e-12,'all')

        error( ...
            'Invalid pressure-volume quadratic discriminant.');

    end


    discriminant = ...
        max(discriminant,0);


    % The "+" root is the physical branch lying between
    % RWC_TLP and 1.

    rwcRaw(idxTurgid) = ...
        ( ...
            -q1 ...
            + sqrt(discriminant) ...
        ) ...
        ./ (2 .* q2);

end


%% ========================================================================
% Physical bounds
% ========================================================================

rwc = ...
    min( ...
        max(rwcRaw,rMin), ...
        1);


lowerBoundActive = ...
    rwcRaw < rMin;


upperBoundActive = ...
    rwcRaw > 1 ...
    | idxSat;


boundActive = ...
    lowerBoundActive ...
    | upperBoundActive;


%% ========================================================================
% Water amounts
% ========================================================================

Wraw = ...
    Wsat ...
    .* rwcRaw;


Wtotal = ...
    Wsat ...
    .* rwc;


%% ------------------------------------------------------------------------
% Apoplastic and symplastic water
%
% af is defined relative to saturated total water content.
% -------------------------------------------------------------------------

Wapo = ...
    af ...
    .* Wsat ...
    .* ones(size(rwc));


Wsym = ...
    Wtotal ...
    - Wapo;


Wsym = ...
    max(Wsym,0);


rwcSym = ...
    ( ...
        rwc ...
        - af ...
    ) ...
    ./ (1 - af);


rwcSym = ...
    min( ...
        max(rwcSym,0), ...
        1);


%% ========================================================================
% Recover osmotic and turgor potentials
% ========================================================================

D = ...
    d0 ...
    + d1 .* rwc;


psiOsm = ...
    A ...
    ./ D;


psiTurgor = ...
    zeros(size(rwc));


idxPressure = ...
    rwc > rTLP;


psiTurgor(idxPressure) = ...
    kTurgor ...
    .* ( ...
        rwc(idxPressure) ...
        - rTLP ...
    );


psiReconstructed = ...
    psiOsm ...
    + psiTurgor;


constitutiveResidual = ...
    psiReconstructed ...
    - psi;


%% ========================================================================
% Effective capacitance
%
% C = dW/dpsi
%   = Wsat / (dpsi/dRWC)
% ========================================================================

dPsiOsm_dRWC = ...
    -A ...
    .* d1 ...
    ./ D.^2;


dPsiTurgor_dRWC = ...
    zeros(size(rwc));


dPsiTurgor_dRWC(idxPressure) = ...
    kTurgor;


dPsi_dRWC = ...
    dPsiOsm_dRWC ...
    + dPsiTurgor_dRWC;


Ceff = ...
    Wsat ...
    ./ dPsi_dRWC;


% No additional storage can be accessed once a numerical
% water-content bound is active.

Ceff(boundActive) = ...
    0;


%% ========================================================================
% Outputs
% ========================================================================

S = ...
    struct;


% Fields expected by the hydraulic-storage interface.

S.amount_mol_m2 = ...
    Wtotal;


S.rawAmount_mol_m2 = ...
    Wraw;


S.capacitanceEffective_mol_m2_MPa = ...
    Ceff;


S.boundActive = ...
    boundActive;


%% P-V diagnostics

S.relativeWaterContent = ...
    rwc;


S.rawRelativeWaterContent = ...
    rwcRaw;


S.relativeWaterContentTurgorLoss = ...
    rTLP;


S.symplasticRelativeWaterContent = ...
    rwcSym;


S.apoplasticFraction = ...
    af;


S.apoplasticWaterAmount_mol_m2 = ...
    Wapo;


S.symplasticWaterAmount_mol_m2 = ...
    Wsym;


S.osmoticPotential_MPa = ...
    psiOsm;


S.turgorPotential_MPa = ...
    psiTurgor;


S.reconstructedPsiLeaf_MPa = ...
    psiReconstructed;


S.constitutiveResidual_MPa = ...
    constitutiveResidual;


S.belowTurgorLoss = ...
    rwc <= rTLP;


S.minimumRelativeWaterContent = ...
    rMin;

end