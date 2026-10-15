function Q = resolve_conditional_prior_mapping(C,M)
%RESOLVE_CONDITIONAL_PRIOR_MAPPING
%
% Determine whether a conditioning context C and mapping metadata M are
% jointly sufficient to permit a numerical conditional prior.
%
% Two conditioning pathways are supported:
%
% A. TARGET-EVIDENCE PATHWAY
%    resolve_conditional_prior_context(C) reports
%    sufficient_for_conditioning for the target dimension.
%
% B. VALIDATED-PREDICTOR PATHWAY
%    the target registry row requires external conditioning AND
%    a validated mapping has complete metadata AND
%    all required predictor fields are present AND
%    the context is inside the declared mapping domain
%    (or extrapolation is allowed).


if nargin<1 || isempty(C)
    C=conditional_prior_context_template();
end
if nargin<2 || isempty(M)
    M=conditional_prior_mapping_template();
end

P=inference_prior_registry();
QC=resolve_conditional_prior_context(C);

Q=struct;
Q.version='2026-09-03_resolve_conditional_prior_mapping_v2';
Q.contextResolution=QC;
Q.mapping=M;

target=string(M.targetDimension);
ids=string({P.id});

Q.targetValid=any(ids==target);

if ~Q.targetValid
    Q.contextEvidenceState="not_applicable";
    Q.targetRequiresExternalConditioning=false;
    Q.requiredPredictorsPresent=false;
    Q.insideDomain=false;
    Q.mappingMetadataComplete=false;
    Q.mappingValidated=false;
    Q.targetEvidencePathway=false;
    Q.validatedPredictorPathway=false;
    Q.sufficiencyPathway="none";
    Q.action="invalid_target_dimension";
    Q.permitted=false;
    return
end

iTarget=find(ids==target,1);
Preg=P(iTarget);

row=QC.dimensionStatus(QC.dimensionStatus.dimension==target,:);
Q.contextEvidenceState=row.evidenceState;
Q.targetRequiresExternalConditioning=logical(Preg.requiresExternalConditioning);

Q.mappingValidated=string(M.validation.status)=="validated";

Q.mappingMetadataComplete = ...
    strlength(string(M.mappingID))>0 && ...
    strlength(string(M.targetParameterization))>0 && ...
    strlength(string(M.mappingMethod))>0 && ...
    strlength(string(M.output.distributionFamily))>0 && ...
    strlength(string(M.output.uncertaintyRepresentation))>0 && ...
    strlength(string(M.calibration.sourceType))>0;

Q.requiredPredictorsPresent = requiredPredictorsPresentLocal( ...
    C,string(M.requiredPredictorFields));

Q.insideDomain = insideDomainLocal(C,M.domain);

Q.targetEvidencePathway = ...
    row.evidenceState=="sufficient_for_conditioning";

Q.validatedPredictorPathway = ...
    Q.targetRequiresExternalConditioning && ...
    Q.mappingMetadataComplete && ...
    Q.mappingValidated && ...
    Q.requiredPredictorsPresent;

Q.sufficiencyPathway="none";

if ~Q.mappingMetadataComplete
    Q.action="mapping_metadata_incomplete";
    Q.permitted=false;

elseif ~Q.mappingValidated
    Q.action="mapping_unvalidated";
    Q.permitted=false;

elseif ~Q.requiredPredictorsPresent
    Q.action="missing_required_predictors";
    Q.permitted=false;

elseif ~Q.insideDomain && ~logical(M.extrapolation.allowed)
    Q.action="outside_mapping_domain";
    Q.permitted=false;

elseif Q.targetEvidencePathway
    Q.sufficiencyPathway="target_evidence";
    Q.action="numerical_mapping_permitted";
    Q.permitted=true;

elseif Q.validatedPredictorPathway
    Q.sufficiencyPathway="validated_predictor_mapping";
    Q.action="numerical_mapping_permitted";
    Q.permitted=true;

else
    Q.action="context_not_sufficient";
    Q.permitted=false;
end
end

function tf=requiredPredictorsPresentLocal(C,fields)

tf=true;

for i=1:numel(fields)
    f=string(fields(i));

    if strlength(f)==0
        continue
    end

    parts=split(f,".");

    if numel(parts)~=2
        error('Predictor field must use block.field form: %s',f);
    end

    block=char(parts(1));
    field=char(parts(2));

    if ~isfield(C,block) || ~isfield(C.(block),field)
        tf=false;
        return
    end

    v=C.(block).(field);

    if isnumeric(v)
        if isempty(v) || ~any(isfinite(v(:)))
            tf=false;
            return
        end
    elseif islogical(v)
        % logical is always a defined value
    else
        if strlength(string(v))==0
            tf=false;
            return
        end
    end
end
end

function tf=insideDomainLocal(C,D)

tf=true;

if ~isempty(D.allowedTaxa)
    vals=string(D.allowedTaxa);
    vals=vals(strlength(vals)>0);
    if ~isempty(vals)
        candidate=[string(C.identity.taxon),string(C.identity.species), ...
            string(C.identity.genus),string(C.identity.family)];
        candidate=candidate(strlength(candidate)>0);
        tf=tf && any(ismember(candidate,vals));
    end
end

if ~isempty(D.allowedPFTs)
    vals=string(D.allowedPFTs);
    vals=vals(strlength(vals)>0);
    if ~isempty(vals)
        tf=tf && ismember(string(C.identity.PFT),vals);
    end
end

if ~isempty(D.allowedGrowthForms)
    vals=string(D.allowedGrowthForms);
    vals=vals(strlength(vals)>0);
    if ~isempty(vals)
        tf=tf && ismember(string(C.identity.growthForm),vals);
    end
end

if ~isempty(D.environmentClasses)
    vals=string(D.environmentClasses);
    vals=vals(strlength(vals)>0);
    if ~isempty(vals)
        tf=tf && ismember(string(C.environment.climateClass),vals);
    end
end
end
