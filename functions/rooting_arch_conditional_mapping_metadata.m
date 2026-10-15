function M = rooting_arch_conditional_mapping_metadata()
%ROOTING_ARCH_CONDITIONAL_MAPPING_METADATA_CANDIDATE
%
% Candidate deterministic constructor for the accepted ROOTING_ARCH
% conditional-mapping lifecycle metadata.


M=conditional_prior_mapping_template();

M.mappingID="ROOTING_ARCH_RPGE_ROBUST_LOGISTIC_NORMAL_V1";
M.targetDimension="ROOTING_ARCH";
M.targetParameterization="three_zone_uptake_fraction_distribution";
M.mappingMethod="robust_huber_ALR_with_tail_calibrated_logistic_normal";
M.mappingVersion="1.0.0";

M.predictorFields=[ ...
"environment.PET_mm"; ...
"environment.MAP_over_PET"; ...
"identity.growthForm"];

M.requiredPredictorFields=M.predictorFields;
M.optionalPredictorFields=strings(0,1);

M.output.distributionFamily="logistic_normal_simplex";
M.output.parameterNames=["zone1_fraction";"zone2_fraction";"zone3_fraction"];
M.output.supportDescription="three positive fractions summing to one";
M.output.uncertaintyRepresentation="2x2 tail-calibrated covariance in ALR space";

M.calibration.sourceType="empirical_compilation";
M.calibration.citation="Schenk & Jackson RPGE / ORNL DAAC rooting profiles";
M.calibration.DOI="10.3334/ORNLDAAC/660";
M.calibration.datasetID="RPGE_D50_D95";
M.calibration.nRecords=547;
M.calibration.taxonomicCoverage="ecosystem profiles spanning multiple vegetation/growth forms";
M.calibration.PFTCoverage="broadleaved trees; needleleaved trees; shrubs; semi-shrubs; grasses; forbs; succulents";
M.calibration.environmentalCoverage="RPGE profiles with PET and MAP available";
M.calibration.notes="M6d-M6h empirical mapping workflow";

M.domain.allowedTaxa=strings(0,1);
M.domain.allowedPFTs=strings(0,1);
M.domain.allowedGrowthForms=[ ...
"broadleavedTrees"; ...
"needleleavedTrees"; ...
"shrubs"; ...
"semiShrubs"; ...
"grasses"; ...
"forbs"; ...
"succulents"];
M.domain.environmentClasses=strings(0,1);
M.domain.notes="No explicit climate-class restriction yet; extrapolation disabled.";

M.validation.status="validated";
M.validation.method="10-fold held-out CV with robust Huber ALR mean and tail-calibrated covariance";
M.validation.checkpoint="historical validation artifact: M6h_RPGE_logistic_normal_uncertainty_calibration_v1";
M.validation.metrics="coverage90=0.934186; coverage95=0.948812; meanZoneRMSE=0.106818; profileP95=0.242403";
M.validation.notes="M6h overall PASS.";

M.extrapolation.allowed=false;
M.extrapolation.policy="";
M.extrapolation.notes="Outside-domain deployment blocked pending explicit validation.";

M.provenance.createdBy="M6i";
M.provenance.createdDate="2026-09-03";
M.provenance.codeVersion="M6i_v3";
M.provenance.notes="Mapping metadata constructed from accepted M6h calibration.";

end
