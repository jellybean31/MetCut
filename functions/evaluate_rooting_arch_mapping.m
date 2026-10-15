function Q = evaluate_rooting_arch_mapping(context,payload)
%EVALUATE_ROOTING_ARCH_MAPPING_PRODUCTION_CANDIDATE
%
% Candidate production implementation for:
%
%   Q = evaluate_rooting_arch_mapping(context,payload)

Q=struct;
Q.version="2026-09-06_rooting_arch_mapping_production_candidate_v1";
Q.status="indeterminate";
Q.zoneUptakeFractions=[NaN NaN NaN];
Q.eta=[NaN NaN];
Q.standardizedPredictors=[];
Q.diagnostics=struct;
Q.diagnostics.simplexClosureResidual=NaN;
Q.diagnostics.finite=false;
Q.diagnostics.reason="";

if ~isstruct(context) || ~isscalar(context)
    Q.diagnostics.reason="invalid_context";
    return
end

if ~isstruct(payload) || ~isscalar(payload)
    Q.diagnostics.reason="invalid_payload";
    return
end

requiredPayload={
    'predictorNames'
    'coefficients'
    'PET_logMean'
    'PET_logSD'
    'aridity_logMean'
    'aridity_logSD'
    'responseParameterization'
    'outputFamily'
};

for i=1:numel(requiredPayload)
    if ~isfield(payload,requiredPayload{i})
        Q.diagnostics.reason="payload_missing_"+string(requiredPayload{i});
        return
    end
end

if string(payload.responseParameterization) ~= ...
        "ALR_f1_over_f3__f2_over_f3"
    Q.diagnostics.reason="unexpected_response_parameterization";
    return
end

if string(payload.outputFamily) ~= "logistic_normal_simplex"
    Q.diagnostics.reason="unexpected_output_family";
    return
end

names=string(payload.predictorNames(:));

expectedNames=[
    "intercept"
    "z_logPET"
    "z_logMAPoverPET"
    "broadleavedTrees"
    "needleleavedTrees"
    "shrubs"
    "semiShrubs"
    "grasses"
    "forbs"
];

if ~isequal(names,expectedNames)
    Q.diagnostics.reason="unexpected_predictor_order";
    return
end

B=double(payload.coefficients);

if ~isequal(size(B),[9 2]) || any(~isfinite(B(:)))
    Q.diagnostics.reason="invalid_coefficient_matrix";
    return
end

if ~isscalar(payload.PET_logMean) || ...
        ~isscalar(payload.PET_logSD) || ...
        ~isfinite(double(payload.PET_logMean)) || ...
        ~isfinite(double(payload.PET_logSD)) || ...
        double(payload.PET_logSD)<=0
    Q.diagnostics.reason="invalid_PET_transform";
    return
end

if ~isscalar(payload.aridity_logMean) || ...
        ~isscalar(payload.aridity_logSD) || ...
        ~isfinite(double(payload.aridity_logMean)) || ...
        ~isfinite(double(payload.aridity_logSD)) || ...
        double(payload.aridity_logSD)<=0
    Q.diagnostics.reason="invalid_aridity_transform";
    return
end

if ~isfield(context,'PET_mm') || ...
        ~isscalar(context.PET_mm) || ...
        ~isfinite(double(context.PET_mm)) || ...
        double(context.PET_mm)<=0
    Q.diagnostics.reason="missing_or_invalid_PET_mm";
    return
end

if ~isfield(context,'MAP_over_PET') || ...
        ~isscalar(context.MAP_over_PET) || ...
        ~isfinite(double(context.MAP_over_PET)) || ...
        double(context.MAP_over_PET)<=0
    Q.diagnostics.reason="missing_or_invalid_MAP_over_PET";
    return
end

x=zeros(1,9);

x(1)=1;
x(2)=(log(double(context.PET_mm))-double(payload.PET_logMean)) / ...
    double(payload.PET_logSD);
x(3)=(log(double(context.MAP_over_PET))-double(payload.aridity_logMean)) / ...
    double(payload.aridity_logSD);

for i=4:9
    fld=char(names(i));

    if ~isfield(context,fld)
        Q.diagnostics.reason="missing_predictor_"+names(i);
        return
    end

    v=context.(fld);

    if ~(isnumeric(v) || islogical(v)) || ...
            ~isscalar(v) || ...
            ~isfinite(double(v))
        Q.diagnostics.reason="invalid_predictor_"+names(i);
        return
    end

    x(i)=double(v);
end

eta=x*B;

Z=[eta 0];
mx=max(Z);
EZ=exp(Z-mx);
f=EZ/sum(EZ);

if any(~isfinite(f)) || any(f<0) || any(f>1)
    Q.diagnostics.reason="invalid_simplex_output";
    return
end

closure=abs(sum(f)-1);

Q.status="resolved";
Q.zoneUptakeFractions=f;
Q.eta=eta;
Q.standardizedPredictors=x;
Q.diagnostics.simplexClosureResidual=closure;
Q.diagnostics.finite=true;
Q.diagnostics.reason="ok";

end
