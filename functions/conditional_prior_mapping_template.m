function M = conditional_prior_mapping_template()
%CONDITIONAL_PRIOR_MAPPING_TEMPLATE
%
% Metadata schema for a numerical conditional-prior mapping.


M=struct;

M.version='2026-09-05_conditional_prior_mapping_template_v2_componentEvidence';

M.mappingID="";
M.targetDimension="";          % e.g. ROOTING_ARCH, LEAF_ANATOMY, I01_I02, W02
M.targetParameterization="";   % explicit output parameterization
M.mappingMethod="";            % e.g. regression, hierarchical model, lookup table
M.mappingVersion="";

M.predictorFields=strings(0,1);
M.requiredPredictorFields=strings(0,1);
M.optionalPredictorFields=strings(0,1);

M.output=struct;
M.output.distributionFamily="";
M.output.parameterNames=strings(0,1);
M.output.supportDescription="";
M.output.uncertaintyRepresentation="";

M.calibration=struct;
M.calibration.sourceType="";
M.calibration.citation="";
M.calibration.DOI="";
M.calibration.datasetID="";
M.calibration.nRecords=NaN;
M.calibration.taxonomicCoverage="";
M.calibration.PFTCoverage="";
M.calibration.environmentalCoverage="";
M.calibration.notes="";

M.domain=struct;
M.domain.allowedTaxa=strings(0,1);
M.domain.allowedPFTs=strings(0,1);
M.domain.allowedGrowthForms=strings(0,1);
M.domain.environmentClasses=strings(0,1);
M.domain.notes="";

M.validation=struct;
M.validation.status="unvalidated"; % allowed: unvalidated / validated / rejected
M.validation.method="";
M.validation.checkpoint="";
M.validation.metrics="";
M.validation.notes="";

M.extrapolation=struct;
M.extrapolation.allowed=false;
M.extrapolation.policy="";
M.extrapolation.notes="";

M.componentEvidence=table; % optional generic component-level evidence/treatment metadata\n\nM.provenance=struct;
M.provenance.createdBy="";
M.provenance.createdDate="";
M.provenance.codeVersion="";
M.provenance.notes="";
end
