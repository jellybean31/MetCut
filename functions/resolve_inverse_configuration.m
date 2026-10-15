function Q = resolve_inverse_configuration(C)
%RESOLVE_INVERSE_CONFIGURATION
%
% Resolve a production inverse configuration into runtime routing.
%
% OUTPUT
% ------
% Q.valid
% Q.inferredDimensions
% Q.marginalizedDimensions
% Q.conditionedDimensions
% Q.fixedDimensions
% Q.disabledOptionalProcesses
% Q.biochemicalMode
% Q.biochemicalLatentDimension
% Q.dimensionRouting
% Q.errors

if nargin<1 || isempty(C)
    C=inverse_configuration_template();
end

P=inference_prior_registry();
ids=string({P.id})';

errors=strings(0,1);

%% Observation validation
obsDesign=string(C.observations.design);
homologues=C.observations.homologues(:)';
obs=C.observations.deltaD_permil(:)';

if strlength(obsDesign)==0
    errors(end+1,1)="observation design missing"; 
end
if isempty(homologues) || isempty(obs) || numel(homologues)~=numel(obs)
    errors(end+1,1)="observation homologues/values invalid"; 
end

%% BIO_MODE12 routing
bioMode="";
bioDim=NaN;

if obsDesign=="C29_only"
    if isequal(homologues,29)
        bioMode="Mode1";
        bioDim=1;
    else
        errors(end+1,1)="C29_only design requires homologue 29 only"; 
    end
elseif obsDesign=="multi_homologue_C27_C29_C31_C33"
    if isequal(homologues,[27 29 31 33])
        bioMode="Mode2";
        bioDim=2;
    else
        errors(end+1,1)="multi-homologue design requires [27 29 31 33]"; 
    end
else
    errors(end+1,1)="unsupported observation design"; 
end

%  user selection must agree with routed mode if provided.
if strlength(string(C.inferenceSelection.biochemicalMode))>0 && ...
        string(C.inferenceSelection.biochemicalMode)~=bioMode
    errors(end+1,1)="biochemical mode conflicts with observation design"; 
end

if isfinite(C.inferenceSelection.biochemicalLatentDimension) && ...
        C.inferenceSelection.biochemicalLatentDimension~=bioDim
    errors(end+1,1)="biochemical latent dimension conflicts with observation design"; 
end

%% Requested inferred dimensions
inferred=string(C.inferenceSelection.dimensions(:));
inferred=inferred(strlength(inferred)>0);

if isempty(inferred) && strlength(bioMode)>0
    inferred="BIO_MODE12";
end

for i=1:numel(inferred)
    j=find(ids==inferred(i),1);
    if isempty(j)
        errors(end+1,1)="unknown inferred dimension: "+inferred(i); 
    elseif ~P(j).freeInferencePermitted
        errors(end+1,1)="free inference not permitted: "+inferred(i); 
    elseif P(j).numericalPriorStatus~="available"
        errors(end+1,1)="available numerical prior absent for inferred dimension: "+inferred(i); 
    end
end

%% Requested marginalized dimensions
marginalized=string(C.marginalization.dimensions(:));
marginalized=marginalized(strlength(marginalized)>0);

for i=1:numel(marginalized)
    j=find(ids==marginalized(i),1);
    if isempty(j)
        errors(end+1,1)="unknown marginalized dimension: "+marginalized(i); 
    elseif ~P(j).marginalizationPermitted
        errors(end+1,1)="marginalization not permitted: "+marginalized(i); 
    elseif P(j).numericalPriorStatus~="available"
        errors(end+1,1)="available numerical prior absent for marginalized dimension: "+marginalized(i); 
    end
end

%% Conditioning context
QC=resolve_conditional_prior_context(C.conditioningContext);

conditioned=strings(0,1);

for target=["ROOTING_ARCH","LEAF_ANATOMY"]
    row=QC.dimensionStatus(QC.dimensionStatus.dimension==target,:);
    if height(row)~=1
        errors(end+1,1)="conditioning resolver missing row: "+target; 
        continue
    end

    if row.evidenceState=="sufficient_for_conditioning"
        conditioned(end+1,1)=target; 
    else
        if target=="ROOTING_ARCH" && ~C.policy.allowMissingROOTING_ARCHContext
            errors(end+1,1)="required ROOTING_ARCH conditioning context missing/insufficient"; 
        elseif target=="LEAF_ANATOMY" && ~C.policy.allowMissingLEAF_ANATOMYContext
            errors(end+1,1)="required LEAF_ANATOMY conditioning context missing/insufficient"; 
        end
    end
end

for target=["I01_I02","W02"]
    row=QC.dimensionStatus(QC.dimensionStatus.dimension==target,:);
    if height(row)==1 && row.evidenceState=="sufficient_for_conditioning"
        conditioned(end+1,1)=target; 
        marginalized(marginalized==target)=[];
    end
end

%% Fixed / optional
fixed=string(C.fixedValues.dimensions(:));
fixed=fixed(strlength(fixed)>0);

disabledOptional=strings(0,1);

optionalIDs=["ROOT_FRAC","SOIL_EVAP_ENRICH","DYNAMIC_RESERVE_C"];
optionalFlags=[ ...
    logical(C.optionalProcesses.ROOT_FRAC), ...
    logical(C.optionalProcesses.SOIL_EVAP_ENRICH), ...
    logical(C.optionalProcesses.DYNAMIC_RESERVE_C)];

for i=1:numel(optionalIDs)
    id=optionalIDs(i);
    if optionalFlags(i)
        errors(end+1,1)="optional process activation not yet authorized in N1 production resolver: "+id; 
    else
        disabledOptional(end+1,1)=id; 
        if ~any(fixed==id)
            fixed(end+1,1)=id; 
        end
    end
end

if C.policy.allowNumericalLEAF_ANATOMYPrior
    errors(end+1,1)="numerical LEAF_ANATOMY prior prohibited"; 
end

allRouted=[inferred; marginalized; conditioned; fixed];
if numel(unique(allRouted))~=numel(allRouted)
    errors(end+1,1)="one or more dimensions assigned to multiple runtime roles"; 
end

runtimeRole=strings(numel(P),1);
status=strings(numel(P),1);

for i=1:numel(P)
    id=P(i).id;

    if any(inferred==id)
        runtimeRole(i)="infer";
        status(i)="active";
    elseif any(marginalized==id)
        runtimeRole(i)="marginalize";
        status(i)="active";
    elseif any(conditioned==id)
        runtimeRole(i)="condition";
        status(i)="active";
    elseif any(fixed==id)
        runtimeRole(i)="fixed";
        status(i)="active";
    elseif P(i).optionalProcess && ~P(i).defaultEnabled
        runtimeRole(i)="disabled_optional";
        status(i)="disabled";
    else
        runtimeRole(i)="UNRESOLVED";
        status(i)="error";
        errors(end+1,1)="dimension unresolved at runtime: "+id; 
    end
end

Q=struct;
Q.version="N1_inverse_configuration_resolver_production_v1";
Q.valid=isempty(errors);
Q.observationDesign=obsDesign;
Q.biochemicalMode=bioMode;
Q.biochemicalLatentDimension=bioDim;
Q.inferredDimensions=unique(inferred,'stable');
Q.marginalizedDimensions=unique(marginalized,'stable');
Q.conditionedDimensions=unique(conditioned,'stable');
Q.fixedDimensions=unique(fixed,'stable');
Q.disabledOptionalProcesses=unique(disabledOptional,'stable');
Q.conditioningResolution=QC;
Q.dimensionRouting=table(ids,runtimeRole,status, ...
    'VariableNames',{'dimension','runtimeRole','status'});
Q.errors=errors;

end
