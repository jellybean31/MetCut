function Q = resolve_conditional_prior_context(C)
%RESOLVE_CONDITIONAL_PRIOR_CONTEXT
%
% Validate external conditioning context against inference_prior_registry().
%
% OUTPUT
% ------
% Q.dimensionStatus
%   One row per registry dimension with:
%     - whether conditioning evidence exists
%     - whether it is direct / borrowed
%     - evidence sufficiency state
%     - permitted action


if nargin<1 || isempty(C)
    C=conditional_prior_context_template();
end

P=inference_prior_registry();

validateContextShapeLocal(C);

rows=cell(numel(P),1);

for i=1:numel(P)

    id=P(i).id;

    hasEvidence=false;
    direct=false;
    borrowed=false;
    state="not_applicable";
    action=P(i).inferenceRole;
    evidenceFields=strings(0,1);

    switch id

        case "ROOTING_ARCH"
            [hasEvidence,direct,borrowed,evidenceFields] = ...
                rootingEvidenceLocal(C);

            if ~hasEvidence
                state="no_evidence";
            elseif direct || numel(evidenceFields)>=2
                state="sufficient_for_conditioning";
            else
                state="partial_evidence";
            end

        case "LEAF_ANATOMY"
            [hasEvidence,direct,borrowed,evidenceFields] = ...
                leafEvidenceLocal(C);

            if ~hasEvidence
                state="no_evidence";
            elseif direct || numel(evidenceFields)>=2
                state="sufficient_for_conditioning";
            else
                state="partial_evidence";
            end

        case "I01_I02"
            [hasEvidence,direct,borrowed,evidenceFields] = ...
                xylemEvidenceLocal(C);

            if ~hasEvidence
                state="no_evidence";
            elseif direct || numel(evidenceFields)>=2
                state="sufficient_for_conditioning";
            else
                state="partial_evidence";
            end

        case "W02"
            [hasEvidence,direct,borrowed,evidenceFields] = ...
                waxEvidenceLocal(C);

            if ~hasEvidence
                state="no_evidence";
            elseif direct || numel(evidenceFields)>=2
                state="sufficient_for_conditioning";
            else
                state="partial_evidence";
            end

        otherwise
            state="not_applicable";
    end

    % Respect explicit borrowing policy.
    if borrowed
        if ~(C.policy.allowTaxonomicBorrowing || ...
             C.policy.allowPFTBorrowing || ...
             C.policy.allowEnvironmentalBorrowing)
            state="partial_evidence";
        end
    end

    if C.policy.requireDirectMeasurement && hasEvidence && ~direct
        state="partial_evidence";
    end

    rows{i}=struct( ...
        'dimension',id, ...
        'registryRole',P(i).inferenceRole, ...
        'priorClass',P(i).priorClass, ...
        'numericalPriorStatus',P(i).numericalPriorStatus, ...
        'requiresExternalConditioning',P(i).requiresExternalConditioning, ...
        'hasEvidence',hasEvidence, ...
        'directEvidence',direct, ...
        'borrowedEvidence',borrowed, ...
        'evidenceState',state, ...
        'evidenceFields',evidenceFields, ...
        'permittedAction',permittedActionLocal(P(i),state));
end

Q=struct;
Q.version='2026-09-03_resolve_conditional_prior_context_v1';
Q.context=C;
Q.registry=P;
Q.dimensionStatus=struct2table(vertcat(rows{:}));

Q.summary=struct;
Q.summary.nDimensions=height(Q.dimensionStatus);
Q.summary.nConditioningReady=nnz( ...
    Q.dimensionStatus.evidenceState=="sufficient_for_conditioning");
Q.summary.nPartial=nnz(Q.dimensionStatus.evidenceState=="partial_evidence");
Q.summary.nNoEvidence=nnz(Q.dimensionStatus.evidenceState=="no_evidence");

Q.designNotes={ ...
    'Evidence sufficiency is a schema/policy decision only.', ...
    'No numerical conditional distribution is produced.', ...
    'Partial evidence must not be silently promoted to a universal prior.', ...
    'Borrowed taxonomic/PFT/environment evidence requires explicit policy permission.', ...
    'Direct measurements can support conditioning but still require a mapping model in a later stage.'};
end

%% =========================================================================
function validateContextShapeLocal(C)

req={'identity','environment','rooting','leafAnatomy','plantSize','wax','provenance','policy'};

for i=1:numel(req)
    if ~isfield(C,req{i})
        error('Conditional context missing required block: %s',req{i});
    end
end
end

%% =========================================================================
function [has,direct,borrowed,fields]=rootingEvidenceLocal(C)

fields=strings(0,1);

if strlength(string(C.rooting.strategyClass))>0
    fields(end+1)="rooting.strategyClass";
end
if isfinite(C.rooting.maxRootingDepth_m)
    fields(end+1)="rooting.maxRootingDepth_m";
end
if isfinite(C.rooting.meanRootingDepth_m)
    fields(end+1)="rooting.meanRootingDepth_m";
end
if isfinite(C.rooting.rootDepth50_m)
    fields(end+1)="rooting.rootDepth50_m";
end
if isfinite(C.rooting.rootDepth95_m)
    fields(end+1)="rooting.rootDepth95_m";
end
if any(isfinite(C.rooting.zoneUptakeFractions))
    fields(end+1)="rooting.zoneUptakeFractions";
end

has=~isempty(fields);
direct=logical(C.rooting.isDirectMeasurement);

borrowed=has && ~direct && any([ ...
    strlength(string(C.identity.PFT))>0, ...
    strlength(string(C.identity.growthForm))>0, ...
    strlength(string(C.provenance.taxonMatch))>0]);
end

%% =========================================================================
function [has,direct,borrowed,fields]=leafEvidenceLocal(C)

fields=strings(0,1);

if isfinite(C.leafAnatomy.leafThickness_um)
    fields(end+1)="leafAnatomy.leafThickness_um";
end
if isfinite(C.leafAnatomy.mesophyllThickness_um)
    fields(end+1)="leafAnatomy.mesophyllThickness_um";
end
if isfinite(C.leafAnatomy.veinDensity_mm_per_mm2)
    fields(end+1)="leafAnatomy.veinDensity_mm_per_mm2";
end
if isfinite(C.leafAnatomy.interveinalDistance_um)
    fields(end+1)="leafAnatomy.interveinalDistance_um";
end
if strlength(string(C.leafAnatomy.cellPackingClass))>0
    fields(end+1)="leafAnatomy.cellPackingClass";
end
if strlength(string(C.leafAnatomy.anatomyClass))>0
    fields(end+1)="leafAnatomy.anatomyClass";
end

has=~isempty(fields);
direct=logical(C.leafAnatomy.isDirectMeasurement);

borrowed=has && ~direct && any([ ...
    strlength(string(C.identity.PFT))>0, ...
    strlength(string(C.identity.taxon))>0, ...
    strlength(string(C.provenance.taxonMatch))>0]);
end

%% =========================================================================
function [has,direct,borrowed,fields]=xylemEvidenceLocal(C)

fields=strings(0,1);

if isfinite(C.plantSize.height_m)
    fields(end+1)="plantSize.height_m";
end
if isfinite(C.plantSize.stemDiameter_cm)
    fields(end+1)="plantSize.stemDiameter_cm";
end
if isfinite(C.plantSize.crownArea_m2)
    fields(end+1)="plantSize.crownArea_m2";
end
if strlength(string(C.plantSize.woodinessClass))>0
    fields(end+1)="plantSize.woodinessClass";
end
if strlength(string(C.identity.growthForm))>0
    fields(end+1)="identity.growthForm";
end

has=~isempty(fields);
direct=false; 

borrowed=has;
end

%% =========================================================================
function [has,direct,borrowed,fields]=waxEvidenceLocal(C)

fields=strings(0,1);

if isfinite(C.wax.leafAge_days)
    fields(end+1)="wax.leafAge_days";
end
if isfinite(C.wax.developmentWindow_days)
    fields(end+1)="wax.developmentWindow_days";
end
if isfinite(C.wax.matureReplacementFraction_perDay)
    fields(end+1)="wax.matureReplacementFraction_perDay";
end
if strlength(string(C.wax.developmentStage))>0
    fields(end+1)="wax.developmentStage";
end

has=~isempty(fields);
direct=logical(C.wax.isDirectMeasurement);
borrowed=has && ~direct;
end

%% =========================================================================
function action=permittedActionLocal(P,state)

switch state
    case "sufficient_for_conditioning"
        if P.requiresExternalConditioning
            action="conditioning_allowed_but_numeric_mapping_required";
        else
            action=P.inferenceRole;
        end

    case "partial_evidence"
        if P.marginalizationPermitted
            action="retain_development_prior_or_marginalize";
        else
            action="no_conditional_prior_yet";
        end

    case "no_evidence"
        if P.marginalizationPermitted
            action="retain_development_prior_or_marginalize";
        elseif P.optionalProcess
            action="retain_default_optional_state";
        else
            action="no_conditional_prior";
        end

    otherwise
        action=P.inferenceRole;
end
end
