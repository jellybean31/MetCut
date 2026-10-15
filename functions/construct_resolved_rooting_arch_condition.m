function Q = construct_resolved_rooting_arch_condition(C,mappingState,payload)
%CONSTRUCT_RESOLVED_ROOTING_ARCH_CONDITION_CANDIDATE
%
% Candidate Helper A implementation.
%
%   Q = construct_resolved_rooting_arch_condition(C,mappingState,payload)

Q=struct;
Q.version="2026-09-06_rooting_helperA_candidate_v2";
Q.status="indeterminate";
Q.zoneUptakeFractions=[NaN NaN NaN];
Q.zoneDepthBounds_m=[0.0 0.3; 0.3 0.6; 0.6 0.9];
Q.sourcePathway="none";
Q.mappingStatus="";
Q.diagnostics=struct;
Q.diagnostics.reason="";
Q.diagnostics.directSimplexClosureResidual=NaN;
Q.diagnostics.mappingAction="";
Q.diagnostics.evaluatorStatus="";

if ~isstruct(C) || ~isscalar(C)
    Q.diagnostics.reason="invalid_C";
    return
end

if ~isstruct(mappingState) || ~isscalar(mappingState)
    Q.diagnostics.reason="invalid_mappingState";
    return
end

%% Direct condition path
if isfield(C,'rooting') && isstruct(C.rooting) && isscalar(C.rooting) && ...
        isfield(C.rooting,'zoneUptakeFractions')

    f=C.rooting.zoneUptakeFractions;

    if isnumeric(f) && numel(f)==3 && any(isfinite(f(:)))
        f=double(reshape(f,1,3));

        if any(~isfinite(f))
            Q.diagnostics.reason="malformed_direct_simplex_nonfinite";
            return
        end

        if any(f<0) || any(f>1)
            Q.diagnostics.reason="malformed_direct_simplex_bounds";
            return
        end

        closure=abs(sum(f)-1);
        Q.diagnostics.directSimplexClosureResidual=closure;

        tol=64*eps(max([1 abs(f)]));

        if closure>tol
            Q.diagnostics.reason="malformed_direct_simplex_closure";
            return
        end

        Q.status="resolved";
        Q.zoneUptakeFractions=f;
        Q.sourcePathway="direct_condition";
        Q.mappingStatus="not_used";
        Q.diagnostics.reason="ok";
        return
    elseif ~isnumeric(f) || numel(f)~=3
        Q.diagnostics.reason="malformed_direct_simplex_shape";
        return
    end
end

%% Numerical mapping permission gate
if ~isfield(mappingState,'action')
    Q.diagnostics.reason="mapping_action_missing";
    return
end

Q.diagnostics.mappingAction=string(mappingState.action);

if string(mappingState.action)~="numerical_mapping_permitted"
    Q.diagnostics.reason="numerical_mapping_not_permitted";
    Q.mappingStatus=string(mappingState.action);
    return
end

if ~isfield(C,'environment') || ...
        ~isstruct(C.environment) || ...
        ~isscalar(C.environment)
    Q.diagnostics.reason="environment_context_missing";
    return
end

if ~isfield(C.environment,'PET_mm') || ...
        ~isscalar(C.environment.PET_mm) || ...
        ~isfinite(double(C.environment.PET_mm)) || ...
        double(C.environment.PET_mm)<=0
    Q.diagnostics.reason="missing_or_invalid_environment_PET_mm";
    return
end

if ~isfield(C.environment,'MAP_over_PET') || ...
        ~isscalar(C.environment.MAP_over_PET) || ...
        ~isfinite(double(C.environment.MAP_over_PET)) || ...
        double(C.environment.MAP_over_PET)<=0
    Q.diagnostics.reason="missing_or_invalid_environment_MAP_over_PET";
    return
end

if ~isfield(C,'identity') || ...
        ~isstruct(C.identity) || ...
        ~isscalar(C.identity) || ...
        ~isfield(C.identity,'growthForm') || ...
        strlength(string(C.identity.growthForm))==0
    Q.diagnostics.reason="missing_growthForm";
    return
end

growthOrder=[
    "broadleavedTrees"
    "needleleavedTrees"
    "shrubs"
    "semiShrubs"
    "grasses"
    "forbs"
];

growthForm=string(C.identity.growthForm);
j=find(growthOrder==growthForm,1);

if isempty(j)
    Q.diagnostics.reason="unsupported_growthForm";
    return
end

ctx=struct;
ctx.PET_mm=double(C.environment.PET_mm);
ctx.MAP_over_PET=double(C.environment.MAP_over_PET);

for i=1:numel(growthOrder)
    ctx.(char(growthOrder(i)))=double(i==j);
end

%% Production numerical evaluator
E=evaluate_rooting_arch_mapping(ctx,payload);

Q.mappingStatus=string(E.status);
Q.diagnostics.evaluatorStatus=string(E.status);

if string(E.status)~="resolved"
    Q.diagnostics.reason="mapping_evaluator_indeterminate";
    return
end

if ~isfield(E,'zoneUptakeFractions') || ...
        ~isnumeric(E.zoneUptakeFractions) || ...
        numel(E.zoneUptakeFractions)~=3
    Q.diagnostics.reason="invalid_mapped_simplex_shape";
    return
end

f=double(reshape(E.zoneUptakeFractions,1,3));

if any(~isfinite(f)) || any(f<0) || any(f>1)
    Q.diagnostics.reason="invalid_mapped_simplex";
    return
end

closure=abs(sum(f)-1);
tol=64*eps(max([1 abs(f)]));

if closure>tol
    Q.diagnostics.reason="mapped_simplex_nonclosure";
    return
end

Q.status="resolved";
Q.zoneUptakeFractions=f;
Q.sourcePathway="validated_mapping";
Q.diagnostics.reason="ok";

end
