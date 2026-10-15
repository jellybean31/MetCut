function [M,B] = rerun_full_model_isotope_source_overrides_with_wax_timing( ...
    Mreference,sourceCfg,stageTraits,accumTraits)
%RERUN_FULL_MODEL_ISOTOPE_SOURCE_OVERRIDES_WITH_WAX_TIMING
%
%
% Extends the validated replay so that the following isotope
% source classes can be perturbed independently while hydraulics remain fixed:
%
%   1. initial soil water by soil zone
%   2. initial xylem lumen storage
%   3. initial xylem living storage
%   4. precipitation isotope forcing by interval
%
% INPUT
% -----
% sourceCfg fields:
%
%   precipDeltaDInterval_permil
%       nTime-1 vector
%
%   initialSoilDeltaDByZone_permil
%       nZone vector
%
%   xylemStorageLumenInitialDeltaD_permil
%       scalar
%
%   xylemStorageLivingInitialDeltaD_permil
%       scalar
%
% All fields are required.
%
% stageTraits / accumTraits:
%   optional.
%
% OUTPUT
% ------
% M includes architecture-v2 wax and standing-wax accumulation.


if nargin<2 || isempty(sourceCfg)
    error('Mreference and sourceCfg are required.');
end

if nargin<3 || isempty(stageTraits)
    stageTraits=wax_er_stage_state_traits();
    stageTraits.enabled=true;
    stageTraits.stage6Delta_permil=0;
    stageTraits.stage7Delta_permil=0;
end

if nargin<4 || isempty(accumTraits)
    accumTraits=wax_accumulation_default_traits();
end

requiredM={ ...
    'timeDays','config','hydraulicCheckpointPath','hydraulics','grid', ...
    'meteorology'};

for k=1:numel(requiredM)
    if ~isfield(Mreference,requiredM{k})
        error('Mreference is missing required field %s.',requiredM{k});
    end
end

requiredS={ ...
    'precipDeltaDInterval_permil', ...
    'initialSoilDeltaDByZone_permil', ...
    'xylemStorageLumenInitialDeltaD_permil', ...
    'xylemStorageLivingInitialDeltaD_permil'};

for k=1:numel(requiredS)
    if ~isfield(sourceCfg,requiredS{k})
        error('sourceCfg is missing required field %s.',requiredS{k});
    end
end

timeDays=Mreference.timeDays(:);
nT=numel(timeDays);
nInt=nT-1;

dP=sourceCfg.precipDeltaDInterval_permil(:);

if numel(dP)~=nInt || any(~isfinite(dP))
    error('sourceCfg.precipDeltaDInterval_permil must contain nTime-1 finite values.');
end

checkpointPath=Mreference.hydraulicCheckpointPath;

if ~isfile(checkpointPath)
    error('Hydraulic checkpoint not found: %s',checkpointPath);
end

Q=load(checkpointPath,'plantCoupled');

if ~isfield(Q,'plantCoupled')
    error('Checkpoint is missing plantCoupled.');
end

plantCoupled=Q.plantCoupled;
Crich=Mreference.hydraulics;
G=Mreference.grid;
met=Mreference.meteorology;
cfg=Mreference.config;

cfg.forcing.precipDeltaD_permil=dP;

zoneIndex=Crich.soil.zoneIndex(:);
nZone=max(zoneIndex);

dZone=sourceCfg.initialSoilDeltaDByZone_permil(:);

if numel(dZone)~=nZone || any(~isfinite(dZone))
    error('initialSoilDeltaDByZone_permil must contain one finite value per soil zone.');
end

cfg.forcing.initialSoilDeltaDByZone_permil=dZone(:).';

deltaDInitialCell_permil=dZone(zoneIndex);

%% Soil isotope transport

Isoil=soil_isotope_transport_richards( ...
    Crich.soil, ...
    deltaDInitialCell_permil, ...
    dP, ...
    [], ...
    cfg.soilIsotope);

%% Root-source boundary

Broot=prepare_root_xylem_isotope_boundary( ...
    Isoil,Crich,struct('Rstandard_DH',cfg.Rstandard_DH));

%% Xylem isotope transport

dL=sourceCfg.xylemStorageLumenInitialDeltaD_permil;
dV=sourceCfg.xylemStorageLivingInitialDeltaD_permil;

if ~isscalar(dL) || ~isfinite(dL) || ~isscalar(dV) || ~isfinite(dV)
    error('Initial xylem storage isotope overrides must be finite scalars.');
end

cfg.xylemIsotope.storageLumenInitialDeltaD=dL;
cfg.xylemIsotope.storageLivingInitialDeltaD=dV;

Xxylem=xylem_transport_dynamic_exactD( ...
    timeDays, ...
    Broot.deltaDNode_permil, ...
    Crich.hydro, ...
    plantCoupled, ...
    cfg.xylemIsotope);

%% Leaf / intracellular / biochemistry

Lleaf=leaf_hydraulics_dynamic_spatial_symplast( ...
    timeDays, ...
    Xxylem.leafBoundary, ...
    Crich.hydro, ...
    met, ...
    cfg.leafIsotope);

Icell=intracellular_water_dynamic_exactD( ...
    timeDays,Lleaf,cfg.intracellular);

Bio=wax_biochemistry_from_intracellular_nominal( ...
    timeDays,Icell,cfg.biochemistry);

%% Assemble trusted replay state

M=struct;
M.version='2026-08-30_H8_full_model_isotope_source_override_v1';
M.timeDays=timeDays;
M.config=cfg;
M.hydraulicCheckpointPath=checkpointPath;
M.hydraulics=Crich;
M.grid=G;
M.meteorology=met;

M.soilIsotopes=Isoil;
M.rootBoundary=Broot;
M.xylemIsotopes=Xxylem;
M.leafIsotopes=Lleaf;
M.intracellular=Icell;
M.biochemistry=Bio;

M.scope=Mreference.scope;
M.knownCurrentOmissions=Mreference.knownCurrentOmissions;

B=build_isotope_fractionation_budget(M,dP);

%% Architecture-v2 + standing wax

M=attach_wax_architecture_v2_to_model(M,stageTraits);

M=attach_wax_accumulation_to_model( ...
    M,accumTraits,'waxStageExpanded');

M.sourceTagging=struct;
M.sourceTagging.version='2026-08-30_H8_source_tagging_v1';
M.sourceTagging.initialSoilDeltaDByZone_permil=dZone;
M.sourceTagging.xylemStorageLumenInitialDeltaD_permil=dL;
M.sourceTagging.xylemStorageLivingInitialDeltaD_permil=dV;
M.sourceTagging.precipDeltaDInterval_permil=dP;
M.sourceTagging.mobileXylemSourceClass= ...
    'inherits initial root source and therefore belongs to initial-soil/root basis';

end
