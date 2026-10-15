function P = forward_observation_adapter(M,C)
%FORWARD_OBSERVATION_ADAPTER
%
% Production forward-observation adapter.
%
% Map an already-computed production model object M into the exact finite
% standing-wax prediction vector required by the likelihood
% schema.
%
% This function:
%   - uses standing wax only;
%   - uses exact model-node time matching only;
%   - performs no interpolation or nearest-neighbor substitution;
%   - blocks empty observation times;
%   - preserves observation-vector order exactly;
%   - fails on nonfinite selected predictions.


if nargin<2
    error('N4c:InputCount','M and C are required.');
end
if ~isstruct(M) || ~isstruct(C)
    error('N4c:InputType','M and C must be structs.');
end
if ~isfield(C,'observations') || ~isstruct(C.observations)
    error('N4c:ObservationSchema','C.observations is required.');
end

O=C.observations;

reqObs={'design','homologues','deltaD_permil','timeDays','isStandingWax'};
for k=1:numel(reqObs)
    if ~isfield(O,reqObs{k})
        error('N4c:ObservationSchema', ...
            'C.observations.%s is required.',reqObs{k});
    end
end

if ~(islogical(O.isStandingWax) || ...
        (isnumeric(O.isStandingWax) && isscalar(O.isStandingWax) && isfinite(O.isStandingWax)))
    error('N4c:StandingWaxFlag','observations.isStandingWax must be scalar logical/numeric.');
end
if ~logical(O.isStandingWax)
    error('N4c:InstantaneousWaxUnsupported', ...
        'Stage N4 currently supports standing-wax inverse observations only.');
end

if ~isfield(M,'biochemistry') || ~isstruct(M.biochemistry) || ...
        ~isfield(M.biochemistry,'waxAccumulation') || ...
        ~isstruct(M.biochemistry.waxAccumulation)
    error('N4c:StandingWaxMissing', ...
        'M.biochemistry.waxAccumulation is required.');
end

A=M.biochemistry.waxAccumulation;
reqA={'poolDeltaD_permil','timeDays','alkaneCarbonNumber'};
for k=1:numel(reqA)
    if ~isfield(A,reqA{k})
        error('N4c:StandingWaxSchema', ...
            'M.biochemistry.waxAccumulation.%s is required.',reqA{k});
    end
end

tModel=A.timeDays(:);
Cmodel=A.alkaneCarbonNumber(:);
Y=A.poolDeltaD_permil;

if isempty(tModel) || any(~isfinite(tModel)) || any(diff(tModel)<=0)
    error('N4c:ModelTime','Standing-wax timeDays must be finite and strictly increasing.');
end
if isempty(Cmodel) || any(~isfinite(Cmodel))
    error('N4c:ModelHomologues','Standing-wax alkaneCarbonNumber must be finite.');
end

nH=numel(Cmodel);
nT=numel(tModel);

if isequal(size(Y),[nH nT])
    %
elseif isequal(size(Y),[nT nH])
    Y=Y.';
else
    error('N4c:StandingWaxDimensions', ...
        'poolDeltaD_permil dimensions must agree with homologues and timeDays.');
end

design=string(O.design);
obs=O.deltaD_permil(:);
hom=O.homologues(:);
tObs=O.timeDays(:);

if isempty(obs)
    error('N4c:ObservationValues','observations.deltaD_permil must be nonempty.');
end
if any(~isfinite(obs))
    error('N4c:ObservationValues','observations.deltaD_permil must be finite.');
end
if isempty(tObs)
    error('N4c:EmptyObservationTime', ...
        'observations.timeDays must be explicit and nonempty in Stage N4.');
end
if any(~isfinite(tObs))
    error('N4c:ObservationTime','observations.timeDays must be finite.');
end

switch design

    case "C29_only"
        expectedHom=29;
        expectedN=1;

    case "multi_homologue_C27_C29_C31_C33"
        expectedHom=[27;29;31;33];
        expectedN=4;

    otherwise
        error('N4c:ObservationDesign', ...
            'Unsupported observation design: %s',design);
end

if numel(obs)~=expectedN
    error('N4c:ObservationLength', ...
        'Observation count does not match observation design.');
end
if ~isequal(hom,expectedHom)
    error('N4c:HomologueOrder', ...
        'Observed homologues do not match the required design/order.');
end
if ~(isscalar(tObs) || numel(tObs)==expectedN)
    error('N4c:ObservationTimeLength', ...
        'timeDays must be scalar or one value per observation.');
end
if isscalar(tObs)
    tPair=repmat(tObs,expectedN,1);
else
    tPair=tObs;
end

hIdx=zeros(expectedN,1);
tIdx=zeros(expectedN,1);
pred=zeros(expectedN,1);

for k=1:expectedN

    ih=find(Cmodel==expectedHom(k));

    if isempty(ih)
        error('N4c:HomologueMissing', ...
            'Required homologue C%d is absent from standing-wax output.',expectedHom(k));
    end
    if numel(ih)~=1
        error('N4c:HomologueDuplicate', ...
            'Required homologue C%d appears more than once.',expectedHom(k));
    end

    it=find(tModel==tPair(k));

    if isempty(it)
        error('N4c:ObservationTimeNotExactNode', ...
            ['Requested observation time %.17g d is not an exact standing-wax ', ...
             'model node. Interpolation and nearest-neighbor mapping are blocked.'], ...
             tPair(k));
    end
    if numel(it)~=1
        error('N4c:ModelTimeDuplicate', ...
            'Standing-wax model time %.17g d is duplicated.',tPair(k));
    end

    hIdx(k)=ih;
    tIdx(k)=it;
    pred(k)=Y(ih,it);
end

if any(~isfinite(pred))
    error('N4c:NonfinitePrediction', ...
        'Selected standing-wax prediction contains nonfinite values.');
end

P=struct;
P.version="N4_forward_observation_adapter_production_v1";
P.predictedDeltaD_permil=pred;
P.observedDeltaD_permil=obs;
P.observationDesign=design;
P.homologues=expectedHom;
P.observationTimeDays=tPair;
P.homologueIndex=hIdx;
P.timeIndex=tIdx;
P.isStandingWax=true;
P.source=struct;
P.source.observableField="M.biochemistry.waxAccumulation.poolDeltaD_permil";
P.source.timeField="M.biochemistry.waxAccumulation.timeDays";
P.source.homologueField="M.biochemistry.waxAccumulation.alkaneCarbonNumber";
P.source.timeMappingPolicy="EXACT_MODEL_NODE_ONLY";
P.source.interpolationUsed=false;
P.source.nearestNeighborUsed=false;
P.source.instantaneousWaxUsed=false;

end
