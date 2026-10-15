function [M,B] = rerun_full_model_temporal_memory_source_case( ...
    Mreference,sourceCfg,memoryCfg,stageTraits,accumTraits)
%RERUN_FULL_MODEL_TEMPORAL_MEMORY_SOURCE_CASE
%
% Replay one fixed-hydraulic temporal-memory/source-tag case.
%
% memoryCfg fields:
%   mobileTransitDays
%   storageExchangeTauDays
%
% The grouped storage exchange timescale is applied to BOTH lumen and living
% storage, matching the established I02 GSA parameterization.
%
% sourceCfg fields:
%   precipDeltaDInterval_permil
%   initialSoilDeltaDByZone_permil
%   xylemStorageLumenInitialDeltaD_permil
%   xylemStorageLivingInitialDeltaD_permil


if nargin<3
    error('Mreference, sourceCfg, and memoryCfg are required.');
end

if nargin<4 || isempty(stageTraits)
    stageTraits=wax_er_stage_state_traits();
    stageTraits.enabled=true;
    stageTraits.stage6Delta_permil=0;
    stageTraits.stage7Delta_permil=0;
end

if nargin<5 || isempty(accumTraits)
    accumTraits=wax_accumulation_default_traits();
end

if ~isfield(memoryCfg,'mobileTransitDays') || ...
        ~isfield(memoryCfg,'storageExchangeTauDays')
    error('memoryCfg requires mobileTransitDays and storageExchangeTauDays.');
end

if ~isscalar(memoryCfg.mobileTransitDays) || ...
        ~isfinite(memoryCfg.mobileTransitDays) || ...
        memoryCfg.mobileTransitDays<=0
    error('mobileTransitDays must be positive and finite.');
end

if ~isscalar(memoryCfg.storageExchangeTauDays) || ...
        ~isfinite(memoryCfg.storageExchangeTauDays) || ...
        memoryCfg.storageExchangeTauDays<=0
    error('storageExchangeTauDays must be positive and finite.');
end

% Start from the validated H8 source-override replay interface, but update
% only the two replayable xylem-memory parameters before transport.

timeDays=Mreference.timeDays(:);
nT=numel(timeDays);
nInt=nT-1;

requiredS={ ...
    'precipDeltaDInterval_permil', ...
    'initialSoilDeltaDByZone_permil', ...
    'xylemStorageLumenInitialDeltaD_permil', ...
    'xylemStorageLivingInitialDeltaD_permil'};

for k=1:numel(requiredS)
    if ~isfield(sourceCfg,requiredS{k})
        error('sourceCfg is missing %s.',requiredS{k});
    end
end

dP=sourceCfg.precipDeltaDInterval_permil(:);
if numel(dP)~=nInt || any(~isfinite(dP))
    error('Invalid precipitation isotope forcing.');
end

checkpointPath=Mreference.hydraulicCheckpointPath;
Q=load(checkpointPath,'plantCoupled');
if ~isfield(Q,'plantCoupled')
    error('Checkpoint is missing plantCoupled.');
end
plantCoupled=Q.plantCoupled;

Crich=Mreference.hydraulics;
met=Mreference.meteorology;
cfg=Mreference.config;

cfg.forcing.precipDeltaD_permil=dP;

zoneIndex=Crich.soil.zoneIndex(:);
nZone=max(zoneIndex);
dZone=sourceCfg.initialSoilDeltaDByZone_permil(:);
if numel(dZone)~=nZone || any(~isfinite(dZone))
    error('Initial soil isotope vector must match soil zones.');
end
cfg.forcing.initialSoilDeltaDByZone_permil=dZone(:).';

deltaDInitialCell_permil=dZone(zoneIndex);

Isoil=soil_isotope_transport_richards( ...
    Crich.soil,deltaDInitialCell_permil,dP,[],cfg.soilIsotope);

Broot=prepare_root_xylem_isotope_boundary( ...
    Isoil,Crich,struct('Rstandard_DH',cfg.Rstandard_DH));

% Apply temporal-memory parameters.
cfg.xylemIsotope.mobileTransitDays=memoryCfg.mobileTransitDays;
cfg.xylemIsotope.storageLumenExchangeTauDays= ...
    memoryCfg.storageExchangeTauDays;
cfg.xylemIsotope.storageLivingExchangeTauDays= ...
    memoryCfg.storageExchangeTauDays;

% Explicit H8 source-tag initial storage isotope states.
cfg.xylemIsotope.storageLumenInitialDeltaD= ...
    sourceCfg.xylemStorageLumenInitialDeltaD_permil;
cfg.xylemIsotope.storageLivingInitialDeltaD= ...
    sourceCfg.xylemStorageLivingInitialDeltaD_permil;

X=xylem_transport_dynamic_exactD( ...
    timeDays,Broot.deltaDNode_permil,Crich.hydro,plantCoupled,cfg.xylemIsotope);

L=leaf_hydraulics_dynamic_spatial_symplast( ...
    timeDays,X.leafBoundary,Crich.hydro,met,cfg.leafIsotope);

I=intracellular_water_dynamic_exactD(timeDays,L,cfg.intracellular);

Bio=wax_biochemistry_from_intracellular_nominal( ...
    timeDays,I,cfg.biochemistry);

M=Mreference;
M.version='2026-08-30_I1_temporal_memory_source_case_v1';
M.config=cfg;
M.soilIsotopes=Isoil;
M.rootBoundary=Broot;
M.xylemIsotopes=X;
M.leafIsotopes=L;
M.intracellular=I;
M.biochemistry=Bio;

B=build_isotope_fractionation_budget(M,dP);

M=attach_wax_architecture_v2_to_model(M,stageTraits);
M=attach_wax_accumulation_to_model(M,accumTraits,'waxStageExpanded');

M.temporalMemoryCase=struct( ...
    'mobileTransitDays',memoryCfg.mobileTransitDays, ...
    'storageExchangeTauDays',memoryCfg.storageExchangeTauDays, ...
    'hydraulicsReused',true);

end
