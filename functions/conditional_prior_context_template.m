function C = conditional_prior_context_template()
%CONDITIONAL_PRIOR_CONTEXT_TEMPLATE
%
% Input schema for external/PFT/taxon/environment conditioning.
%
% This function creates EMPTY CONTEXT METADATA.
%
% The schema is designed for the inference_prior_registry() dimensions:
%   ROOTING_ARCH
%   LEAF_ANATOMY
%   I01_I02
%   W02

C=struct;

C.version='2026-09-03_conditional_prior_context_template_v1';

%% ========================================================================
% Observation / specimen identity
% ========================================================================

C.identity=struct;
C.identity.contextID="";
C.identity.taxon="";
C.identity.species="";
C.identity.genus="";
C.identity.family="";
C.identity.PFT="";
C.identity.growthForm="";
C.identity.lifeHistory="";
C.identity.siteID="";
C.identity.sampleID="";

%% ========================================================================
% Environment
% ========================================================================

C.environment=struct;
C.environment.climateClass="";
C.environment.meanAnnualPrecip_mm=NaN;
C.environment.meanAnnualTemperature_C=NaN;
C.environment.aridityIndex=NaN;
C.environment.soilProfileClass="";
C.environment.soilDepth_m=NaN;
C.environment.waterTableDepth_m=NaN;
C.environment.eventState="";
C.environment.antecedentMoistureState="";
C.environment.notes="";

%% ========================================================================
% Rooting information
% ========================================================================

C.rooting=struct;
C.rooting.strategyClass="";
C.rooting.maxRootingDepth_m=NaN;
C.rooting.meanRootingDepth_m=NaN;
C.rooting.rootDepth50_m=NaN;
C.rooting.rootDepth95_m=NaN;
C.rooting.zoneUptakeFractions=[NaN NaN NaN];
C.rooting.zoneDepthBounds_m=[NaN NaN; NaN NaN; NaN NaN];
C.rooting.isDirectMeasurement=false;
C.rooting.measurementMethod="";
C.rooting.temporalContext="";
C.rooting.environmentSpecific=false;
C.rooting.notes="";

%% ========================================================================
% Leaf anatomy information
% ========================================================================

C.leafAnatomy=struct;
C.leafAnatomy.leafThickness_um=NaN;
C.leafAnatomy.mesophyllThickness_um=NaN;
C.leafAnatomy.veinDensity_mm_per_mm2=NaN;
C.leafAnatomy.interveinalDistance_um=NaN;
C.leafAnatomy.bundleSheathFraction=NaN;
C.leafAnatomy.cellPackingClass="";
C.leafAnatomy.anatomyClass="";
C.leafAnatomy.isDirectMeasurement=false;
C.leafAnatomy.measurementMethod="";
C.leafAnatomy.notes="";

%% ========================================================================
% Plant size / xylem conditioning
% ========================================================================

C.plantSize=struct;
C.plantSize.height_m=NaN;
C.plantSize.stemDiameter_cm=NaN;
C.plantSize.crownArea_m2=NaN;
C.plantSize.woodinessClass="";
C.plantSize.ageClass="";
C.plantSize.notes="";

%% ========================================================================
% Wax development / turnover
% ========================================================================

C.wax=struct;
C.wax.leafAge_days=NaN;
C.wax.developmentStage="";
C.wax.developmentWindow_days=NaN;
C.wax.matureReplacementFraction_perDay=NaN;
C.wax.isDirectMeasurement=false;
C.wax.notes="";

%% ========================================================================
% Provenance / evidence
% ========================================================================

C.provenance=struct;
C.provenance.sourceType="";          % e.g. direct_measurement, literature,
                                     % database, expert_assignment
C.provenance.citation="";
C.provenance.DOI="";
C.provenance.datasetID="";
C.provenance.recordID="";
C.provenance.year=NaN;
C.provenance.evidenceQuality="";     % high / medium / low / unknown
C.provenance.geographicMatch="";     % exact / regional / biome / none
C.provenance.taxonMatch="";          % species / genus / family / PFT / none
C.provenance.environmentMatch="";    % exact / similar / generic / none
C.provenance.notes="";

%% ========================================================================
% User-specified conditioning policy
% ========================================================================

C.policy=struct;
C.policy.allowTaxonomicBorrowing=false;
C.policy.allowPFTBorrowing=false;
C.policy.allowEnvironmentalBorrowing=false;
C.policy.requireDirectMeasurement=false;
C.policy.minimumEvidenceQuality="";
C.policy.notes="";
end
