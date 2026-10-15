function Ustate = gsa_build_upstream_state(Mreference,sample)
%GSA_BUILD_UPSTREAM_STATE Rerun only xylem -> leaf -> intracellular stages.
if nargin<2, sample=struct; end
S=gsa_apply_isotope_sample(Mreference,sample);
t=Mreference.timeDays(:); Crich=Mreference.hydraulics; met=Mreference.meteorology;
Broot=Mreference.rootBoundary; cfg=S.cfg;
Q=load(Mreference.hydraulicCheckpointPath,'plantCoupled');
if ~isfield(Q,'plantCoupled'), error('Checkpoint is missing plantCoupled.'); end
plantCoupled=Q.plantCoupled;
cfg.xylemIsotope.storageLumenInitialDeltaD=Broot.initialDeltaD_permil;
cfg.xylemIsotope.storageLivingInitialDeltaD=Broot.initialDeltaD_permil;
X=xylem_transport_dynamic_exactD(t,Broot.deltaDNode_permil,Crich.hydro,plantCoupled,cfg.xylemIsotope);
L=leaf_hydraulics_dynamic_spatial_symplast(t,X.leafBoundary,Crich.hydro,met,cfg.leafIsotope);
I=intracellular_water_dynamic_exactD(t,L,cfg.intracellular);
Mup=Mreference; Mup.version='2026-08-28_GSA_upstream_state_v1'; Mup.config=cfg;
Mup.xylemIsotopes=X; Mup.leafIsotopes=L; Mup.intracellular=I; Mup.gsaUpstreamOverrides=S;
Ustate=struct('Mup',Mup,'xylem',X,'leaf',L,'intracellular',I,'cfg',cfg);
end
