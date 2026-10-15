function [M,B,O] = run_isotope_gsa_case(Mreference,sample)
%RUN_ISOTOPE_GSA_CASE Replay one non-hydraulic GSA case.
%
% Hydraulics, soil isotope transport, and root-source mixing are reused
% from Mreference. The function reruns xylem -> leaf -> intracellular ->
% biochemistry with the requested overrides. 
%
% Sampled parameters are defined by
% gsa_isotope_parameter_definitions.

if nargin<2, sample=struct; end
S=gsa_apply_isotope_sample(Mreference,sample);

t=Mreference.timeDays(:);
Crich=Mreference.hydraulics;
met=Mreference.meteorology;
Broot=Mreference.rootBoundary;
cfg=S.cfg;

checkpointPath=Mreference.hydraulicCheckpointPath;
Q=load(checkpointPath,'plantCoupled');
if ~isfield(Q,'plantCoupled')
    error('Checkpoint is missing plantCoupled.');
end
plantCoupled=Q.plantCoupled;

% Xylem replay from the unchanged root-source boundary.
cfg.xylemIsotope.storageLumenInitialDeltaD=Broot.initialDeltaD_permil;
cfg.xylemIsotope.storageLivingInitialDeltaD=Broot.initialDeltaD_permil;
X=xylem_transport_dynamic_exactD( ...
    t,Broot.deltaDNode_permil,Crich.hydro,plantCoupled,cfg.xylemIsotope);

% Spatial leaf replay.
L=leaf_hydraulics_dynamic_spatial_symplast( ...
    t,X.leafBoundary,Crich.hydro,met,cfg.leafIsotope);

% Intracellular replay.
I=intracellular_water_dynamic_exactD(t,L,cfg.intracellular);

% Parameterized biochemical bridge.
Bio=wax_biochemistry_from_intracellular_parameters( ...
    t,I,cfg.biochemistry,S.sourceVals,S.reactionAlpha);

M=Mreference;
M.version='2026-08-28_isotope_GSA_replay_v1';
M.config=cfg;
M.xylemIsotopes=X;
M.leafIsotopes=L;
M.intracellular=I;
M.biochemistry=Bio;
M.gsaOverrides=S;

% Preserve original precipitation forcing when rebuilding the budget.
dP=Mreference.config.forcing.precipDeltaD_permil;
if isscalar(dP)
    dP=repmat(dP,numel(t)-1,1);
else
    dP=dP(:);
end
B=build_isotope_fractionation_budget(M,dP);
O=summarize_gsa_case(M,B);
end
