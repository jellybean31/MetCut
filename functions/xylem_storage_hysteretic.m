function S = xylem_storage_hysteretic( ...
    psiX, psiPrevious, lumenPrevious, traits, hasPrevious)
%XYLEM_STORAGE_HYSTERETIC
%
% Two-component xylem storage with:
%
%   1. hysteretic lumen/capillary storage
%   2. reversible living/elastic storage
%
%
% Lumen storage:
%
%   During drying:
%       can only stay constant or decrease.
%       Release follows the release pressure curve.
%
%   During rewetting:
%       can only stay constant or increase.
%       Refill follows a separate refill pressure curve.
%
%
% Living storage:
%
%   Remains reversible and pressure-dependent.
%
%
% hasPrevious = false:
%
%   Used at model initialization. Lumen storage is initialized
%   according to the drying/release equilibrium curve.


%% ========================================================================
% Required traits
% ========================================================================

required = { ...
    'storageLumenWaterInitial_mol_m2', ...
    'livingStorageWaterInitial_mol_m2', ...
    'storageLumenReleaseStart_MPa', ...
    'storageLumenReleaseEnd_MPa', ...
    'storageLumenRefillStart_MPa', ...
    'storageLumenRefillFull_MPa', ...
    'livingStoragePsiFull_MPa', ...
    'livingStoragePsiEmpty_MPa'};


for k = 1:numel(required)

    if ~isfield(traits,required{k})

        error( ...
            'traits.%s is required for hysteretic storage.', ...
            required{k});

    end

end

releaseActive = ...
    false;

refillActive = ...
    false;

%% ========================================================================
% Initial capacities
% ========================================================================

S0Lumen = ...
    traits.storageLumenWaterInitial_mol_m2;


S0Living = ...
    traits.livingStorageWaterInitial_mol_m2;


%% ========================================================================
% Lumen pressure ranges
% ========================================================================

psiReleaseStart = ...
    traits.storageLumenReleaseStart_MPa;


psiReleaseEnd = ...
    traits.storageLumenReleaseEnd_MPa;


psiRefillStart = ...
    traits.storageLumenRefillStart_MPa;


psiRefillFull = ...
    traits.storageLumenRefillFull_MPa;


%% ========================================================================
% Living-storage pressure range
% ========================================================================

psiLivingFull = ...
    traits.livingStoragePsiFull_MPa;


psiLivingEmpty = ...
    traits.livingStoragePsiEmpty_MPa;


%% ========================================================================
% Validation
% ========================================================================

if psiReleaseStart <= psiReleaseEnd

    error([ ...
        'storageLumenReleaseStart_MPa must be less negative than ' ...
        'storageLumenReleaseEnd_MPa.']);

end


if psiRefillFull <= psiRefillStart

    error([ ...
        'storageLumenRefillFull_MPa must be less negative than ' ...
        'storageLumenRefillStart_MPa.']);

end


if psiRefillStart <= psiReleaseStart

    error([ ...
        'For hysteresis, storageLumenRefillStart_MPa should be ' ...
        'less negative than storageLumenReleaseStart_MPa.']);

end


if psiLivingFull <= psiLivingEmpty

    error([ ...
        'livingStoragePsiFull_MPa must be less negative than ' ...
        'livingStoragePsiEmpty_MPa.']);

end


%% ========================================================================
% Living / elastic storage
%
% Reversible for now.
% ========================================================================

fLiving = ...
    linearFraction( ...
        psiX, ...
        psiLivingFull, ...
        psiLivingEmpty);


SLiving = ...
    S0Living ...
    .* fLiving;


CLiving = ...
    0;


if ...
        psiX < psiLivingFull ...
        && psiX > psiLivingEmpty

    CLiving = ...
        S0Living ...
        ./ (psiLivingFull - psiLivingEmpty);

end


%% ========================================================================
% Lumen / capillary storage
% ========================================================================

releaseFraction = ...
    linearFraction( ...
        psiX, ...
        psiReleaseStart, ...
        psiReleaseEnd);


releaseTarget = ...
    S0Lumen ...
    .* releaseFraction;


refillFraction = ...
    linearRefillFraction( ...
        psiX, ...
        psiRefillStart, ...
        psiRefillFull);


refillTarget = ...
    S0Lumen ...
    .* refillFraction;


%% ------------------------------------------------------------------------
% Initialization
% -------------------------------------------------------------------------

if ...
        ~hasPrevious ...
        || ~isfinite(psiPrevious) ...
        || ~isfinite(lumenPrevious)

    SLumen = ...
        releaseTarget;


    if ...
            psiX < psiReleaseStart ...
            && psiX > psiReleaseEnd

        CLumen = ...
            S0Lumen ...
            ./ (psiReleaseStart - psiReleaseEnd);

    else

        CLumen = ...
            0;

    end


    branchCode = ...
        0;


%% ------------------------------------------------------------------------
% Drying
%
% psi becomes more negative.
%
% Storage is not allowed to increase.
% -------------------------------------------------------------------------

elseif psiX < psiPrevious

    SLumen = ...
        min( ...
            lumenPrevious, ...
            releaseTarget);


    releaseActive = ...
           releaseTarget < lumenPrevious ...
        && psiX < psiReleaseStart ...
        && psiX > psiReleaseEnd;


    if releaseActive

        CLumen = ...
            S0Lumen ...
            ./ (psiReleaseStart - psiReleaseEnd);

    else

        CLumen = ...
            0;

    end


    branchCode = ...
        1;


%% ------------------------------------------------------------------------
% Rewetting
%
% psi becomes less negative.
%
% Storage is not allowed to decrease.
% -------------------------------------------------------------------------

elseif psiX > psiPrevious

    SLumen = ...
        max( ...
            lumenPrevious, ...
            refillTarget);


    refillActive = ...
           refillTarget > lumenPrevious ...
        && psiX > psiRefillStart ...
        && psiX < psiRefillFull;


    if refillActive

        CLumen = ...
            S0Lumen ...
            ./ (psiRefillFull - psiRefillStart);

    else

        CLumen = ...
            0;

    end


    branchCode = ...
        -1;


%% ------------------------------------------------------------------------
% No pressure change
% -------------------------------------------------------------------------

else

    SLumen = ...
        lumenPrevious;


    CLumen = ...
        0;


    branchCode = ...
        0;

end


%% ========================================================================
% Numerical bounds
% ========================================================================

SLumen = ...
    min( ...
        max(SLumen,0), ...
        S0Lumen);


SLiving = ...
    min( ...
        max(SLiving,0), ...
        S0Living);


%% ========================================================================
% Outputs
% ========================================================================

S.total_mol_m2 = ...
    SLumen ...
    + SLiving;


S.lumen_mol_m2 = ...
    SLumen;


S.living_mol_m2 = ...
    SLiving;


S.lumenFractionRemaining = ...
    SLumen ...
    ./ S0Lumen;


S.livingFractionRemaining = ...
    SLiving ...
    ./ S0Living;


S.totalFractionRemaining = ...
    S.total_mol_m2 ...
    ./ (S0Lumen + S0Living);


S.capacitanceLumen_mol_m2_MPa = ...
    CLumen;


S.capacitanceLiving_mol_m2_MPa = ...
    CLiving;


S.capacitanceTotal_mol_m2_MPa = ...
    CLumen ...
    + CLiving;


% Branch code:
%
%    1 = drying/release
%   -1 = wetting/refill
%    0 = initialization/no pressure change

S.lumenBranchCode = ...
    branchCode;


S.releaseTarget_mol_m2 = ...
    releaseTarget;


S.refillTarget_mol_m2 = ...
    refillTarget;

S.lumenReleaseActive = ...
    releaseActive;

S.lumenRefillActive = ...
    refillActive;

end



%% =========================================================================
function f = linearFraction( ...
    psi, psiFull, psiEmpty)
%LINEARFRACTION
%
% 1 at psi >= psiFull
% 0 at psi <= psiEmpty

f = ...
    (psi - psiEmpty) ...
    ./ (psiFull - psiEmpty);


f = ...
    min( ...
        max(f,0), ...
        1);

end



%% =========================================================================
function f = linearRefillFraction( ...
    psi, psiStart, psiFull)
%LINEARREFILLFRACTION
%
% 0 at psi <= psiStart
% 1 at psi >= psiFull

f = ...
    (psi - psiStart) ...
    ./ (psiFull - psiStart);


f = ...
    min( ...
        max(f,0), ...
        1);

end