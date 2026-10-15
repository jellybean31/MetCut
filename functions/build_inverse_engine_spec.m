function S = build_inverse_engine_spec(C)
%BUILD_INVERSE_ENGINE_SPEC
%
% Production engine-state specification builder.
%
% Resolves one production inverse configuration and converts registry
% inference roles into sampler-vector packing.
%
% OUTPUT
% ------
% S.Q
% S.biochemical
% S.nuisance
% S.stochasticVectorLabels
% S.dimension


if nargin~=1 || ~isstruct(C) || ~isscalar(C)
    error('N6:InvalidConfiguration','C must be one scalar struct.');
end

Q=resolve_inverse_configuration(C);
if ~Q.valid
    msg=strjoin(string(Q.errors(:)),' | ');
    error('N6:InvalidResolvedConfiguration', ...
        'resolve_inverse_configuration returned invalid: %s',msg);
end

P=inference_prior_registry();
if isstruct(P), P=struct2table(P); end
if ~istable(P)
    error('N6:InvalidRegistry','inference_prior_registry must return table/struct.');
end

runtimeIDs=["B08";"W01";"I01";"I02";"I07";"I08";"W02"];
groupIDs=["B08_W01";"I01_I02";"I07_I08";"W02"];

% Expand registry groups exactly from parameterIDs.
groupMap=struct;
for i=1:numel(groupIDs)
    j=find(string(P.id)==groupIDs(i),1);
    if isempty(j)
        error('N6:MissingRegistryGroup','Missing registry group %s.',groupIDs(i));
    end
    groupMap.(char(groupIDs(i)))=normalizeIDsLocal(P.parameterIDs(j));
end

inferred=string(Q.inferredDimensions(:));
marginalized=string(Q.marginalizedDimensions(:));
conditioned=string(Q.conditionedDimensions(:));
fixed=string(Q.fixedDimensions(:));

% BIO routing.
bioRole="";
if any(inferred=="BIO_MODE12")
    bioRole="INFERRED";
elseif any(marginalized=="BIO_MODE12")
    bioRole="MARGINALIZED";
elseif any(fixed=="BIO_MODE12")
    error('N6:FixedBiochemicalUnsupported', ...
        'Fixed BIO_MODE12 numerical payload has not been contracted.');
elseif any(conditioned=="BIO_MODE12")
    error('N6:ConditionedBiochemicalUnsupported', ...
        'Conditioned BIO_MODE12 is not supported.');
else
    error('N6:BiochemicalUnresolved','BIO_MODE12 has no runtime owner.');
end

bioDim=double(Q.biochemicalLatentDimension);
if ~(isscalar(bioDim) && isfinite(bioDim) && any(bioDim==[1 2]))
    error('N6:InvalidBiochemicalDimension','BIO latent dimension must be 1 or 2.');
end

% Nuisance role expansion.
role=strings(numel(runtimeIDs),1);
sourceGroup=strings(numel(runtimeIDs),1);

for ig=1:numel(groupIDs)
    gid=groupIDs(ig);
    ids=groupMap.(char(gid));

    if any(inferred==gid)
        r="INFERRED";
    elseif any(marginalized==gid)
        r="MARGINALIZED";
    elseif any(fixed==gid)
        r="FIXED";
    elseif any(conditioned==gid)
        r="CONDITIONED";
    else
        r="";
    end

    for k=1:numel(ids)
        j=find(runtimeIDs==ids(k),1);
        if isempty(j)
            error('N6:UnexpectedRuntimeID', ...
                'Registry group %s contains unexpected runtime ID %s.',gid,ids(k));
        end
        if strlength(role(j))>0
            error('N6:DuplicateOwnership','Runtime ID %s has duplicate ownership.',ids(k));
        end
        role(j)=r;
        sourceGroup(j)=gid;
    end
end

if any(strlength(role)==0)
    error('N6:MissingOwnership', ...
        'Every runtime nuisance coordinate must have exactly one owner.');
end

for gid=["I01_I02","W02"]
    if any(conditioned==gid)
        error('N6:ConditionedNuisanceMappingUnavailable', ...
            ['%s is routed to CONDITIONED, but no production numerical ' ...
             'conditioning map to runtime nuisance coordinates is yet available.'],gid);
    end
end

% Fixed physical values must exist per individual runtime ID.
fixedPhysical=struct;
for i=1:numel(runtimeIDs)
    if role(i)=="FIXED"
        id=runtimeIDs(i);
        if ~isfield(C,'fixedValues') || ~isfield(C.fixedValues,'values') || ...
                ~isstruct(C.fixedValues.values) || ...
                ~isfield(C.fixedValues.values,char(id))
            error('N6:MissingFixedValue', ...
                'C.fixedValues.values.%s is required in physical space.',id);
        end
        x=C.fixedValues.values.(char(id));
        if ~isnumeric(x) || ~isscalar(x) || ~isfinite(x)
            error('N6:InvalidFixedValue', ...
                'C.fixedValues.values.%s must be one finite numeric scalar.',id);
        end
        T=apply_nuisance_transform(id,x,"physical_to_sampler");
        if ~isscalar(T.valid) || ~logical(T.valid) || ~isscalar(T.value) || ~isfinite(T.value)
            error('N6:FixedValueOutsideSupport', ...
                'Fixed physical value for %s is outside support.',id);
        end
        fixedPhysical.(char(id))=double(x);
    end
end

stochasticNuisance=runtimeIDs(role=="INFERRED" | role=="MARGINALIZED");

labels=strings(0,1);
for k=1:bioDim
    labels(end+1,1)="BIO_z"+k; 
end
labels=[labels; stochasticNuisance];

S=struct;
S.version="N6_engine_spec_production_v1";
S.Q=Q;
S.registry=P;
S.groupMap=groupMap;
S.biochemical=struct( ...
    'role',bioRole, ...
    'dimension',bioDim, ...
    'vectorIndices',(1:bioDim));
S.nuisance=table(runtimeIDs,sourceGroup,role, ...
    'VariableNames',{'parameterID','registryGroup','role'});
S.fixedPhysicalValues=fixedPhysical;
S.stochasticNuisanceIDs=stochasticNuisance;
S.stochasticVectorLabels=labels;
S.dimension=numel(labels);
S.nuisanceVectorIndices=struct;

offset=bioDim;
for i=1:numel(stochasticNuisance)
    S.nuisanceVectorIndices.(char(stochasticNuisance(i)))=offset+i;
end

end

%% =========================================================================
function ids=normalizeIDsLocal(raw)

if iscell(raw), raw=raw{1}; end
if isstring(raw)
    ids=raw(:).';
elseif ischar(raw)
    ids=string({raw});
elseif iscellstr(raw)
    ids=string(raw(:)).';
else
    ids=string(raw(:)).';
end
ids=ids(strlength(ids)>0);

end
