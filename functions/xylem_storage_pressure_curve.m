function S = xylem_storage_pressure_curve( ...
    psiX, traits)
%XYLEM_STORAGE_PRESSURE_CURVE
%
% Pressure-dependent two-component xylem water storage.
%
% Components:
%
%   1. lumen / capillary storage
%   2. living / elastic storage
%
% Each component is represented initially by a reversible piecewise
% linear storage-pressure relation:
%
%                S0                    psi >= psiFull
%
%   S(psi) =     linear                psiEmpty < psi < psiFull
%
%                0                     psi <= psiEmpty
%



%% ========================================================================
% Required fields
% ========================================================================

required = { ...
    'storageLumenWaterInitial_mol_m2', ...
    'livingStorageWaterInitial_mol_m2', ...
    'storageLumenPsiFull_MPa', ...
    'storageLumenPsiEmpty_MPa', ...
    'livingStoragePsiFull_MPa', ...
    'livingStoragePsiEmpty_MPa'};


for k = 1:numel(required)

    if ~isfield(traits,required{k})

        error( ...
            'traits.%s is required for two-component storage.', ...
            required{k});

    end

end


%% ========================================================================
% Initial inventories
% ========================================================================

S0Lumen = ...
    traits.storageLumenWaterInitial_mol_m2;


S0Living = ...
    traits.livingStorageWaterInitial_mol_m2;


%% ========================================================================
% Pressure ranges
% ========================================================================

psiFullLumen = ...
    traits.storageLumenPsiFull_MPa;


psiEmptyLumen = ...
    traits.storageLumenPsiEmpty_MPa;


psiFullLiving = ...
    traits.livingStoragePsiFull_MPa;


psiEmptyLiving = ...
    traits.livingStoragePsiEmpty_MPa;


%% ========================================================================
% Validate
%
% Example:
%
%   psiFull  = -0.1 MPa
%   psiEmpty = -1.5 MPa
%
% Therefore:
%
%   psiFull > psiEmpty
% ========================================================================

if psiFullLumen <= psiEmptyLumen

    error([ ...
        'storageLumenPsiFull_MPa must be less negative than ' ...
        'storageLumenPsiEmpty_MPa.']);

end


if psiFullLiving <= psiEmptyLiving

    error([ ...
        'livingStoragePsiFull_MPa must be less negative than ' ...
        'livingStoragePsiEmpty_MPa.']);

end


%% ========================================================================
% Piecewise storage fractions
% ========================================================================

fLumen = ...
    storageFractionLinear( ...
        psiX, ...
        psiFullLumen, ...
        psiEmptyLumen);


fLiving = ...
    storageFractionLinear( ...
        psiX, ...
        psiFullLiving, ...
        psiEmptyLiving);


%% ========================================================================
% Component amounts
% ========================================================================

SLumen = ...
    S0Lumen ...
    .* fLumen;


SLiving = ...
    S0Living ...
    .* fLiving;


%% ========================================================================
% Effective capacitance
%
% dS/dPsi is positive:
%
% storage increases as water potential becomes less negative.
% ========================================================================

CLumen = ...
    zeros(size(psiX));


activeLumen = ...
       psiX < psiFullLumen ...
    &  psiX > psiEmptyLumen;


CLumen(activeLumen) = ...
    S0Lumen ...
    ./ (psiFullLumen - psiEmptyLumen);


CLiving = ...
    zeros(size(psiX));


activeLiving = ...
       psiX < psiFullLiving ...
    &  psiX > psiEmptyLiving;


CLiving(activeLiving) = ...
    S0Living ...
    ./ (psiFullLiving - psiEmptyLiving);


%% ========================================================================
% Output
% ========================================================================

S.total_mol_m2 = ...
    SLumen + SLiving;


S.lumen_mol_m2 = ...
    SLumen;


S.living_mol_m2 = ...
    SLiving;


S.lumenFractionRemaining = ...
    fLumen;


S.livingFractionRemaining = ...
    fLiving;


S.totalFractionRemaining = ...
    S.total_mol_m2 ...
    ./ (S0Lumen + S0Living);


S.capacitanceLumen_mol_m2_MPa = ...
    CLumen;


S.capacitanceLiving_mol_m2_MPa = ...
    CLiving;


S.capacitanceTotal_mol_m2_MPa = ...
    CLumen + CLiving;

end



%% =========================================================================
function f = storageFractionLinear( ...
    psi, psiFull, psiEmpty)
%STORAGEFRACTIONLINEAR

f = ...
    (psi - psiEmpty) ...
    ./ (psiFull - psiEmpty);


f = ...
    min( ...
        max(f,0), ...
        1);

end