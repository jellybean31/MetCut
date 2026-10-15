function C = inverse_configuration_template()
%INVERSE_CONFIGURATION_TEMPLATE
%
% Production inverse configuration schema.

C=struct;

C.version="N1_inverse_configuration_production_v1";

C.observations=struct;
C.observations.design="";
C.observations.homologues=[];
C.observations.deltaD_permil=[];
C.observations.timeDays=[];
C.observations.isStandingWax=true;

C.observationError=struct;
C.observationError.sigma_permil=[];
C.observationError.covariance=[];
C.observationError.errorModel="independent_gaussian";
C.observationError.notes="";

C.conditioningContext=conditional_prior_context_template();

C.inferenceSelection=struct;
C.inferenceSelection.dimensions=strings(0,1);
C.inferenceSelection.biochemicalMode="";
C.inferenceSelection.biochemicalLatentDimension=NaN;
C.inferenceSelection.notes="";

C.marginalization=struct;
C.marginalization.dimensions=strings(0,1);
C.marginalization.method="unspecified";
C.marginalization.drawCount=NaN;
C.marginalization.notes="";

C.fixedValues=struct;
C.fixedValues.dimensions=strings(0,1);
C.fixedValues.values=struct;
C.fixedValues.notes="";

C.optionalProcesses=struct;
C.optionalProcesses.ROOT_FRAC=false;
C.optionalProcesses.SOIL_EVAP_ENRICH=false;
C.optionalProcesses.DYNAMIC_RESERVE_C=false;
C.optionalProcesses.notes="";

C.modelDiscrepancy=struct;
C.modelDiscrepancy.enabled=true;
C.modelDiscrepancy.observationSD_permil=[];
C.modelDiscrepancy.surrogateCovariance=[];
C.modelDiscrepancy.modelCovariance=[];
C.modelDiscrepancy.notes="";

C.forwardAdapter=struct;
C.forwardAdapter.hydraulicCheckpoint="";
C.forwardAdapter.referenceModelCheckpoint="";
C.forwardAdapter.forwardMode="production_full_model";
C.forwardAdapter.useCachedHydraulics=true;
C.forwardAdapter.notes="";

C.diagnostics=struct;
C.diagnostics.requireFiniteLogPosterior=true;
C.diagnostics.requireSupportCompliance=true;
C.diagnostics.requireConditionalContext=true;
C.diagnostics.requireForwardSuccess=true;
C.diagnostics.notes="";

C.policy=struct;
C.policy.registryFreezeID="StageM_InferenceArchitecture_v1";
C.policy.allowUnsupportedFreeInference=false;
C.policy.allowUnsupportedMarginalization=false;
C.policy.allowNumericalLEAF_ANATOMYPrior=false;
C.policy.allowMissingROOTING_ARCHContext=false;
C.policy.allowMissingLEAF_ANATOMYContext=false;

C.provenance=struct;
C.provenance.createdBy="Stage N production";
C.provenance.createdDate="2026-09-05";
C.provenance.notes="Production inverse configuration schema.";
end
