function Q = resolve_rooting_arch_hydraulic_compatibility( ...
    resolvedRootingCondition,hydraulicReferenceSummary,hydraulicReferenceIdentity)
%RESOLVE_ROOTING_ARCH_HYDRAULIC_COMPATIBILITY
%
% Production ROOTING_ARCH hydraulic compatibility resolver.
%
% INPUTS
% ------
% resolvedRootingCondition
%   struct with:
%     zoneUptakeFractions : 1x3 structural root-fraction simplex
%     zoneDepthBounds_m   : 3x2 zone bounds [0 .3; .3 .6; .6 .9]
%
% hydraulicReferenceSummary
%   struct with:
%     threeZoneStructuralRootFractions : 1x3 structural root-fraction simplex
%     zoneDepthBounds_m                 : 3x2 accepted zone bounds
%
% hydraulicReferenceIdentity
%   arbitrary provenance value preserved in output.
%
% OUTPUT
% ------
% Q.status = compatible | incompatible | indeterminate


acceptedBounds=[0 0.3;0.3 0.6;0.6 0.9];

Q=struct;
Q.version="2026-09-06_rooting_arch_hydraulic_compatibility_v1";
Q.status="indeterminate";
Q.compatible=false;
Q.reason="";

Q.conditionFractions=nan(1,3);
Q.referenceFractions=nan(1,3);
Q.differenceVector=nan(1,3);

Q.comparisonMetric="max_abs_fraction_difference";
Q.metricValue=NaN;
Q.toleranceValue=0;
Q.toleranceSource="PRE-O13 exact/machine-precision identity policy";

Q.conditionFieldsUsed=["zoneUptakeFractions","zoneDepthBounds_m"];
Q.hydraulicReferenceFieldsUsed= ...
    ["threeZoneStructuralRootFractions","zoneDepthBounds_m"];

Q.hydraulicReferenceIdentity=hydraulicReferenceIdentity;
Q.mappingVersion="three_zone_structural_root_fraction_simplex_v1";

Q.diagnostics=struct;
Q.diagnostics.conditionPresent=false;
Q.diagnostics.referencePresent=false;
Q.diagnostics.conditionGeometryValid=false;
Q.diagnostics.referenceGeometryValid=false;
Q.diagnostics.conditionSimplexValid=false;
Q.diagnostics.referenceSimplexValid=false;
Q.diagnostics.machinePrecisionTolerance=NaN;
Q.diagnostics.proceedCachedHydraulics=false;
Q.diagnostics.requiresHydraulicRebuildOrAlternateReference=false;
Q.diagnostics.failureClass="configuration_error";
Q.diagnostics.posteriorMinusInf=false;
Q.diagnostics.supportRejection=false;
Q.diagnostics.automaticFallback=false;
Q.diagnostics.automaticRebuild=false;
Q.diagnostics.automaticReferenceSelection=false;
Q.diagnostics.mutationPerformed=false;

if ~isstruct(resolvedRootingCondition) || ...
        ~isfield(resolvedRootingCondition,'zoneUptakeFractions') || ...
        ~isfield(resolvedRootingCondition,'zoneDepthBounds_m')
    Q.reason="condition_representation_missing";
    return
end
Q.diagnostics.conditionPresent=true;

if ~isstruct(hydraulicReferenceSummary) || ...
        ~isfield(hydraulicReferenceSummary,'threeZoneStructuralRootFractions') || ...
        ~isfield(hydraulicReferenceSummary,'zoneDepthBounds_m')
    Q.reason="hydraulic_reference_representation_missing";
    return
end
Q.diagnostics.referencePresent=true;

conditionBounds=resolvedRootingCondition.zoneDepthBounds_m;
referenceBounds=hydraulicReferenceSummary.zoneDepthBounds_m;

Q.diagnostics.conditionGeometryValid= ...
    validBoundsLocal(conditionBounds,acceptedBounds);

if ~Q.diagnostics.conditionGeometryValid
    Q.reason="condition_zone_geometry_missing_or_not_definition_equivalent";
    return
end

Q.diagnostics.referenceGeometryValid= ...
    validBoundsLocal(referenceBounds,acceptedBounds);

if ~Q.diagnostics.referenceGeometryValid
    Q.reason="reference_zone_geometry_missing_or_not_definition_equivalent";
    return
end

conditionFractions=resolvedRootingCondition.zoneUptakeFractions;
referenceFractions=hydraulicReferenceSummary.threeZoneStructuralRootFractions;

Q.diagnostics.conditionSimplexValid=validSimplexLocal(conditionFractions);

if ~Q.diagnostics.conditionSimplexValid
    Q.reason="condition_simplex_missing_or_invalid";
    return
end

Q.diagnostics.referenceSimplexValid=validSimplexLocal(referenceFractions);

if ~Q.diagnostics.referenceSimplexValid
    Q.reason="reference_simplex_missing_or_invalid";
    return
end

conditionFractions=reshape(conditionFractions,1,3);
referenceFractions=reshape(referenceFractions,1,3);

Q.conditionFractions=conditionFractions;
Q.referenceFractions=referenceFractions;
Q.differenceVector=conditionFractions-referenceFractions;
Q.metricValue=max(abs(Q.differenceVector));

if all(conditionFractions==referenceFractions)
    Q.status="compatible";
    Q.compatible=true;
    Q.reason="exact_structural_fraction_identity";
    Q.diagnostics.proceedCachedHydraulics=true;
    return
end

scale=max([1 abs(conditionFractions) abs(referenceFractions)]);
machineTol=8*eps(max(scale));
Q.diagnostics.machinePrecisionTolerance=machineTol;

if Q.metricValue<=machineTol
    Q.status="compatible";
    Q.compatible=true;
    Q.reason="difference_consistent_with_machine_precision_roundoff";
    Q.diagnostics.proceedCachedHydraulics=true;
    return
end

Q.status="incompatible";
Q.compatible=false;
Q.reason="nonzero_structural_fraction_difference_without_validated_tolerance";
Q.diagnostics.proceedCachedHydraulics=false;
Q.diagnostics.requiresHydraulicRebuildOrAlternateReference=true;

end

function tf=validSimplexLocal(x)
tf=false;
if ~isnumeric(x) || numel(x)~=3 || any(~isfinite(x))
    return
end
x=reshape(x,1,3);
if any(x<0) || any(x>1)
    return
end
tol=16*eps(max(1,max(abs(x))));
tf=abs(sum(x)-1)<=tol;
end

function tf=validBoundsLocal(x,accepted)
tf=false;
if ~isnumeric(x) || ~isequal(size(x),[3 2]) || any(~isfinite(x(:)))
    return
end
scale=max([1;abs(x(:));abs(accepted(:))]);
tol=8*eps(max(scale));
tf=max(abs(x(:)-accepted(:)))<=tol;
end
