function E = inverse_runtime_dispatcher(Mreference,C,proposal)
%INVERSE_RUNTIME_DISPATCHER
%
% Production runtime dispatcher with LEAF_ANATOMY application.
%
% INPUT
% -----
% Mreference
%   Validated reference model with hydraulics.
%
% C
%   Production inverse configuration.
%
% proposal
%   Scalar struct with:
%
%     .biochemicalLatent
%         Mode1: one finite scalar z1.
%         Mode2: two finite values [z1 z2].
%
%     .nuisanceSampler
%         Scalar struct containing sampler-space values for:
%           B08, W01, I01, I02, I07, I08, W02
%
% OUTPUT
% ------
% E
%   Runtime evaluation record. On valid support:
%     .prediction
%     .logLikelihood
%     .likelihoodDetail
%     .model
%     .budget
%     .resolvedConfiguration
%     .biochemicalTransform
%     .nuisance
%     .physicalSample
%     .stageTraits
%     .accumulationTraits
%

if nargin~=3
    error('N4i:InputCount', ...
        'Expected Mreference, C, and proposal.');
end

if ~isstruct(Mreference) || ~isscalar(Mreference)
    error('N4i:InvalidReferenceModel', ...
        'Mreference must be a scalar struct.');
end

if ~isstruct(C) || ~isscalar(C)
    error('N4i:InvalidConfiguration', ...
        'C must be a scalar struct.');
end

if ~isstruct(proposal) || ~isscalar(proposal)
    error('N4i:InvalidProposal', ...
        'proposal must be a scalar struct.');
end

%% ------------------------------------------------------------------------
% Resolve production inverse configuration
% -------------------------------------------------------------------------

Q=resolve_inverse_configuration(C);

if ~isfield(Q,'valid') || ~logical(Q.valid)
    error('N4i:InvalidResolvedConfiguration', ...
        'Production inverse configuration did not resolve as valid.');
end

if ~isfield(C,'forwardAdapter') || ...
        ~isfield(C.forwardAdapter,'useCachedHydraulics') || ...
        ~isscalar(C.forwardAdapter.useCachedHydraulics) || ...
        ~logical(C.forwardAdapter.useCachedHydraulics)
    error('N4i:CachedHydraulicsRequired', ...
        'N4i requires forwardAdapter.useCachedHydraulics=true.');
end

if ~isfield(C.forwardAdapter,'forwardMode') || ...
        string(C.forwardAdapter.forwardMode)~="production_full_model"
    error('N4i:UnsupportedForwardMode', ...
        'N4i supports forwardMode=production_full_model only.');
end

if ~isfield(C,'observations') || ...
        ~isfield(C.observations,'isStandingWax') || ...
        ~isscalar(C.observations.isStandingWax) || ...
        ~logical(C.observations.isStandingWax)
    error('N4i:StandingWaxRequired', ...
        'N4i supports standing-wax observations only.');
end

%% ------------------------------------------------------------------------
% ROOTING_ARCH structural compatibility gate
% -------------------------------------------------------------------------
rootingMetadata=rooting_arch_conditional_mapping_metadata();
rootingPayload=rooting_arch_mapping_payload();

rootingMappingState=resolve_conditional_prior_mapping( ...
    C.conditioningContext,rootingMetadata);


rootingCondition=construct_resolved_rooting_arch_condition( ...
    C.conditioningContext,rootingMappingState,rootingPayload);

if ~isfield(rootingCondition,'status') || ...
        string(rootingCondition.status)~="resolved"
    error('N4i:ROOTINGARCHConditionUnresolved', ...
        'ROOTING_ARCH condition could not be resolved before replay.');
end


rootingReference=construct_rooting_arch_hydraulic_reference_summary(Mreference);

if ~isfield(rootingReference,'status') || ...
        string(rootingReference.status)~="resolved"
    error('N4i:ROOTINGARCHHydraulicReferenceUnresolved', ...
        'Cached-hydraulic ROOTING_ARCH structural reference could not be resolved.');
end


rootingCompatibilityReference=struct;
rootingCompatibilityReference.threeZoneStructuralRootFractions= ...
    rootingReference.zoneUptakeFractions;
rootingCompatibilityReference.zoneDepthBounds_m= ...
    rootingReference.zoneDepthBounds_m;

rootingCompatibility=resolve_rooting_arch_hydraulic_compatibility( ...
    rootingCondition,rootingCompatibilityReference, ...
    rootingReference.hydraulicCheckpointPath);

if ~isfield(rootingCompatibility,'status')
    error('N4i:ROOTINGARCHCompatibilityUnresolved', ...
        'ROOTING_ARCH hydraulic compatibility status was not returned.');
end

rootingCompatibilityStatus=string(rootingCompatibility.status);

if rootingCompatibilityStatus=="incompatible"
    error('N4i:ROOTINGARCHHydraulicIncompatible', ...
        ['Conditioned ROOTING_ARCH is structurally incompatible with the ' ...
         'cached hydraulic reference. Hydraulic rebuild/reparameterization ' ...
         'is not permitted in this runtime.']);
elseif rootingCompatibilityStatus~="compatible"
    error('N4i:ROOTINGARCHCompatibilityUnresolved', ...
        'ROOTING_ARCH hydraulic compatibility is indeterminate.');
end

%% ------------------------------------------------------------------------
% Initialize output
% -------------------------------------------------------------------------

E=struct;
E.version="N4_runtime_dispatcher_production_v1";
E.supportValid=false;
E.rejectionReason="";
E.logLikelihood=-Inf;
E.prediction=[];
E.likelihoodDetail=[];
E.model=[];
E.budget=[];
E.resolvedConfiguration=Q;
E.biochemicalTransform=[];
E.nuisance=[];
E.physicalSample=[];
E.stageTraits=[];
E.accumulationTraits=[];
E.leafAnatomyApplication=[];
E.logPriorBiochemicalUnnormalized=-Inf;
E.logAbsNuisanceJacobian=NaN;

%% ------------------------------------------------------------------------
% BIO_MODE12 routing and transform
% -------------------------------------------------------------------------

if ~isfield(proposal,'biochemicalLatent')
    error('N4i:MissingBiochemicalLatent', ...
        'proposal.biochemicalLatent is required.');
end

z=double(proposal.biochemicalLatent(:).');

if any(~isfinite(z))
    E.rejectionReason="NONFINITE_BIOCHEMICAL_LATENT";
    return
end

mode=string(Q.biochemicalMode);
dim=double(Q.biochemicalLatentDimension);

switch mode
    case "Mode1"
        if dim~=1 || numel(z)~=1
            error('N4i:Mode1DimensionMismatch', ...
                'Mode1 requires exactly one biochemical latent coordinate.');
        end

        % N4i0b accepted exact embedding.
        zRuntime=[z(1) 0];

    case "Mode2"
        if dim~=2 || numel(z)~=2
            error('N4i:Mode2DimensionMismatch', ...
                'Mode2 requires exactly two biochemical latent coordinates.');
        end

        zRuntime=z;

    otherwise
        error('N4i:UnsupportedBiochemicalMode', ...
            'Unsupported biochemical mode: %s',mode);
end

Tbio=biochemical_latent_mode12_transform(zRuntime);

E.biochemicalTransform=Tbio;
E.biochemicalLatentInput=z;
E.biochemicalLatentRuntime=zRuntime;
E.logPriorBiochemicalUnnormalized=Tbio.logPriorUnnormalized(1);

if ~isscalar(Tbio.feasible) || ~logical(Tbio.feasible)
    E.rejectionReason="BIO_MODE12_OUTSIDE_B01_B07_SUPPORT";
    return
end

B=Tbio.B(1,:);

if numel(B)~=7 || any(~isfinite(B))
    error('N4i:InvalidBiochemicalTransformOutput', ...
        'BIO_MODE12 transform did not return one finite B01-B07 vector.');
end

%% ------------------------------------------------------------------------
% Nuisance transforms
% -------------------------------------------------------------------------

if ~isfield(proposal,'nuisanceSampler') || ...
        ~isstruct(proposal.nuisanceSampler) || ...
        ~isscalar(proposal.nuisanceSampler)
    error('N4i:MissingNuisanceSampler', ...
        'proposal.nuisanceSampler must be a scalar struct.');
end

nuisanceIDs=["B08","W01","I01","I02","I07","I08","W02"];
physical=struct;
transformRecord=struct;
logJac=0;

for k=1:numel(nuisanceIDs)
    id=nuisanceIDs(k);

    if ~isfield(proposal.nuisanceSampler,char(id))
        error('N4i:MissingNuisanceCoordinate', ...
            'proposal.nuisanceSampler.%s is required.',id);
    end

    value=proposal.nuisanceSampler.(char(id));

    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
        E.rejectionReason="INVALID_NUISANCE_SAMPLER_"+id;
        return
    end

    Tk=apply_nuisance_transform(id,value,"sampler_to_physical");
    transformRecord.(char(id))=Tk;

    if ~isscalar(Tk.valid) || ~logical(Tk.valid) || ...
            ~isscalar(Tk.value) || ~isfinite(Tk.value)
        E.rejectionReason="NUISANCE_OUTSIDE_SUPPORT_"+id;
        E.nuisance=transformRecord;
        return
    end

    physical.(char(id))=Tk.value;

    if ~isscalar(Tk.logAbsDxDz) || ~isfinite(Tk.logAbsDxDz)
        error('N4i:InvalidNuisanceJacobian', ...
            'Nuisance transform %s returned an invalid Jacobian.',id);
    end

    logJac=logJac+Tk.logAbsDxDz;
end

E.nuisance=transformRecord;
E.nuisancePhysical=physical;
E.logAbsNuisanceJacobian=logJac;

%% ------------------------------------------------------------------------
% Assemble one source-backed physical sample
% -------------------------------------------------------------------------

sample=struct;

for k=1:7
    sample.(sprintf('B%02d',k))=B(k);
end

sample.B08=physical.B08;
sample.I01=physical.I01;
sample.I02=physical.I02;
sample.I07=physical.I07;
sample.I08=physical.I08;

E.physicalSample=sample;

%% ------------------------------------------------------------------------
% Cached-hydraulic isotope replay
% -------------------------------------------------------------------------

[M,Bbudget]=run_isotope_gsa_case(Mreference,sample);

%% ------------------------------------------------------------------------
% LEAF_ANATOMY component-wise runtime conditioning
% -------------------------------------------------------------------------

[M,Bbudget,~,leafAnatomyApplication]= ...
    apply_leaf_anatomy_condition( ...
        M,C.conditioningContext);

E.leafAnatomyApplication=leafAnatomyApplication;

%% ------------------------------------------------------------------------
% Architecture-v2 biochemical stage states
% -------------------------------------------------------------------------

stageTraits=wax_er_stage_state_traits();
stageTraits.enabled=true;
stageTraits.stage6Delta_permil=0;
stageTraits.stage7Delta_permil=0;

M=attach_wax_architecture_v2_to_model(M,stageTraits);

E.stageTraits=stageTraits;

%% ------------------------------------------------------------------------
% Standing-wax accumulation
% -------------------------------------------------------------------------

accumTraits=wax_accumulation_default_traits();
accumTraits.productionMode="developmentPlusMaintenance";
accumTraits.developmentWindowDays=physical.W01;
accumTraits.matureReplacementFractionPerDay=physical.W02;
accumTraits.lossRatePerDay=physical.W02;

M=attach_wax_accumulation_to_model( ...
    M,accumTraits,'waxStageExpanded');

E.accumulationTraits=accumTraits;

%% ------------------------------------------------------------------------
% Exact configured observation extraction
% -------------------------------------------------------------------------

P=forward_observation_adapter(M,C);

if ~isfield(P,'predictedDeltaD_permil') || ...
        isempty(P.predictedDeltaD_permil) || ...
        any(~isfinite(P.predictedDeltaD_permil),'all')
    error('N4i:InvalidForwardPrediction', ...
        'forward_observation_adapter returned an invalid prediction.');
end


L=resolve_observation_likelihood_schema( ...
    C,P.predictedDeltaD_permil);

[ll,detail]=gaussian_observation_likelihood( ...
    L.observedDeltaD_permil, ...
    L.predictedDeltaD_permil, ...
    L.sigma_permil, ...
    [], ...
    L.discrepancyCovariance);

if ~isscalar(ll) || ~isfinite(ll)
    error('N4i:InvalidLikelihood', ...
        'Gaussian observation likelihood returned a nonfinite value.');
end

%% ------------------------------------------------------------------------
% Finalize
% -------------------------------------------------------------------------

E.supportValid=true;
E.rejectionReason="";
E.logLikelihood=ll;
E.prediction=P;
E.likelihoodSchema=L;
E.likelihoodDetail=detail;
E.model=M;
E.budget=Bbudget;

end

