function [M,B,O,A] = apply_leaf_anatomy_condition(Mbase,context)
%APPLY_LEAF_ANATOMY_CONDITION_PREO49_CANDIDATE
% Production LEAF_ANATOMY runtime application.
%
% Applies model-geometry evidence:
%   leafAnatomy.interveinalDistance_um
%   leafAnatomy.mesophyllThickness_um
%
% Reuses authoritative hydraulics and xylem unchanged, then recomputes:
%   leaf -> intracellular -> biochemistry -> budget/summary.


if nargin~=2
    error('N4i:LEAFANATOMYInputCount','Expected Mbase context.');
end
if ~isstruct(Mbase) || ~isscalar(Mbase)
    error('N4i:LEAFANATOMYInvalidBaseModel','Mbase must be a scalar struct.');
end
if ~isstruct(context) || ~isscalar(context)
    error('N4i:LEAFANATOMYInvalidContext','context must be a scalar struct.');
end

requiredModelFields={'timeDays','hydraulics','meteorology','config','xylemIsotopes','biochemistry'};
for i=1:numel(requiredModelFields)
    if ~isfield(Mbase,requiredModelFields{i})
        error('N4i:LEAFANATOMYMissingBaseField','Mbase.%s is required.',requiredModelFields{i});
    end
end
if ~isfield(Mbase.config,'leafIsotope') || ~isfield(Mbase.config.leafIsotope,'anatomy')
    error('N4i:LEAFANATOMYMissingReferenceAnatomy','Reference leaf anatomy is required.');
end
if ~isfield(Mbase.config,'intracellular') || ~isfield(Mbase.config,'biochemistry')
    error('N4i:LEAFANATOMYMissingDownstreamConfig','Intracellular/biochemistry config required.');
end
if ~isfield(context,'leafAnatomy') || ~isstruct(context.leafAnatomy)
    error('N4i:LEAFANATOMYLeafAnatomyContextMissing','context.leafAnatomy is required.');
end

Q=resolve_conditional_prior_context(context);
row=Q.dimensionStatus(Q.dimensionStatus.dimension=="LEAF_ANATOMY",:);
if height(row)~=1
    error('N4i:LEAFANATOMYLeafStatusRowCount','Expected exactly one LEAF_ANATOMY row.');
end
evidenceState=string(row.evidenceState);
if evidenceState=="no_evidence"
    error('N4i:LEAFANATOMYConditionUnresolved','LEAF_ANATOMY evidence is absent.');
elseif evidenceState=="partial_evidence"
    error('N4i:LEAFANATOMYConditionPartial','LEAF_ANATOMY evidence is partial.');
elseif evidenceState~="sufficient_for_conditioning"
    error('N4i:LEAFANATOMYConditionStateUnsupported', ...
        'Unsupported LEAF_ANATOMY evidence state: %s',evidenceState);
end

LA=context.leafAnatomy;
forbiddenFields={'effectiveFlowAreaFraction','effectiveFlowAreaFraction_value', ...
    'qLeaf','qLeaf_m_inv','tortuosity','lateralPathFraction', ...
    'verticalPathFraction','flowFraction'};
for i=1:numel(forbiddenFields)
    f=forbiddenFields{i};
    if isfield(LA,f) && ~isempty(LA.(f))
        error('N4i:LEAFANATOMYUnauthorizedLeafAnatomyField', ...
            'Numerical LEAF_ANATOMY conditioning of %s is not authorized.',f);
    end
end

cfg=Mbase.config;
refAnatomy=cfg.leafIsotope.anatomy;
newAnatomy=refAnatomy;
appliedFields=strings(0,1);
appliedValues=[];


if isfield(LA,'interveinalDistance_um') && ...
        isnumeric(LA.interveinalDistance_um) && ...
        isscalar(LA.interveinalDistance_um) && ...
        isfinite(LA.interveinalDistance_um)

    if LA.interveinalDistance_um<=0
        error('N4i:LEAFANATOMYInvalidInterveinalDistance', ...
            'interveinalDistance_um must be positive when numerically supplied.');
    end
    v=double(LA.interveinalDistance_um)*1e-6;
    newAnatomy.interveinalDistance_m=v;
    appliedFields(end+1,1)="interveinalDistance_m"; 
    appliedValues(end+1,1)=v; 
end

if isfield(LA,'mesophyllThickness_um') && ...
        isnumeric(LA.mesophyllThickness_um) && ...
        isscalar(LA.mesophyllThickness_um) && ...
        isfinite(LA.mesophyllThickness_um)

    if LA.mesophyllThickness_um<=0
        error('N4i:LEAFANATOMYInvalidMesophyllThickness', ...
            'mesophyllThickness_um must be positive when numerically supplied.');
    end
    v=double(LA.mesophyllThickness_um)*1e-6;
    newAnatomy.mesophyllThickness_m=v;
    appliedFields(end+1,1)="mesophyllThickness_m"; 
    appliedValues(end+1,1)=v; 
end

if isempty(appliedFields)
    % Compatibility / forward-reference branch:
    % sufficient conditioning evidence does not necessarily imply that the
    % runtime layer has an authorized translated numerical geometry value.
    % In that case retain the complete forward-reference anatomy unchanged,
    % perform no raw-trait translation, and continue without a leaf replay.
    M=Mbase;

    dP=Mbase.config.forcing.precipDeltaD_permil;
    t=Mbase.timeDays(:);
    if isscalar(dP)
        dP=repmat(dP,numel(t)-1,1);
    else
        dP=dP(:);
    end
    B=build_isotope_fractionation_budget(M,dP);
    O=summarize_gsa_case(M,B);

    A=struct;
    A.version='PREO53a_leaf_anatomy_compatibility_candidate_v1';
    A.evidenceState=evidenceState;
    A.applicationStatus="reference_retained_no_applicable_geometry";
    A.geometryApplied=false;
    A.appliedFields=appliedFields;
    A.appliedValues_SI=appliedValues;
    A.referenceAnatomy=refAnatomy;
    A.conditionedAnatomy=refAnatomy;
    A.hydraulicsReused=true;
    A.xylemReused=true;
    A.leafReplayPerformed=false;
    A.numericalQLeafPriorConstructed=false;
    A.aEffInferenceOrMarginalizationUsed=false;
    A.rawTraitTranslationPerformed=false;
    A.policy="sufficient evidence with no authorized numerical geometry retains forward-reference anatomy";
    return
end

cfg.leafIsotope.anatomy=newAnatomy;
t=Mbase.timeDays(:);
Crich=Mbase.hydraulics;
met=Mbase.meteorology;
X=Mbase.xylemIsotopes;

L=leaf_hydraulics_dynamic_spatial_symplast( ...
    t,X.leafBoundary,Crich.hydro,met,cfg.leafIsotope);
I=intracellular_water_dynamic_exactD(t,L,cfg.intracellular);

if ~isfield(Mbase.biochemistry,'sourceParameterValues') || ...
        ~isfield(Mbase.biochemistry,'reactionAlphaValues')
    error('N4i:LEAFANATOMYBiochemicalParameterStateMissing', ...
        'Mbase biochemistry parameter values are required for exact downstream replay.');
end
Bio=wax_biochemistry_from_intracellular_parameters( ...
    t,I,cfg.biochemistry,Mbase.biochemistry.sourceParameterValues, ...
    Mbase.biochemistry.reactionAlphaValues);

M=Mbase;
M.version='2026-09-07_PREO51_leaf_anatomy_condition_candidate_v1';
M.config=cfg;
M.leafIsotopes=L;
M.intracellular=I;
M.biochemistry=Bio;

dP=Mbase.config.forcing.precipDeltaD_permil;
if isscalar(dP)
    dP=repmat(dP,numel(t)-1,1);
else
    dP=dP(:);
end
B=build_isotope_fractionation_budget(M,dP);
O=summarize_gsa_case(M,B);

A=struct;
A.version='PREO53a_leaf_anatomy_compatibility_candidate_v1';
A.evidenceState=evidenceState;
A.applicationStatus="explicit_geometry_applied";
A.geometryApplied=true;
A.appliedFields=appliedFields;
A.appliedValues_SI=appliedValues;
A.referenceAnatomy=refAnatomy;
A.conditionedAnatomy=newAnatomy;
A.hydraulicsReused=true;
A.xylemReused=true;
A.leafReplayPerformed=true;
A.numericalQLeafPriorConstructed=false;
A.aEffInferenceOrMarginalizationUsed=false;
A.rawTraitTranslationPerformed=false;
A.policy="component-wise empirical geometry application with structural forward-reference remainder";
end

function tf=isfiniteScalarLocal(x)
tf=isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x);
end
