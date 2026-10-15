function L = resolve_observation_likelihood_schema(C,prediction)
%RESOLVE_OBSERVATION_LIKELIHOOD_SCHEMA
%
% Production adapter from the inverse configuration schema
% to Gaussian observation likelihood.
%
% Production interface. 
%
% Output L contains:
%   observedDeltaD_permil
%   predictedDeltaD_permil
%   sigma_permil
%   discrepancyCovariance
%   errorModel
%   observationDesign
%   homologues
%   timeDays
%   isStandingWax
%   sourcePolicy

if nargin<2
    error('Configuration and prediction are required.');
end

if ~isstruct(C)
    error('C must be a configuration struct.');
end

reqTop={'observations','observationError','modelDiscrepancy'};
for k=1:numel(reqTop)
    if ~isfield(C,reqTop{k})
        error('Configuration missing required block %s.',reqTop{k});
    end
end

O=C.observations;
E=C.observationError;
D=C.modelDiscrepancy;

reqObs={'design','homologues','deltaD_permil','timeDays','isStandingWax'};
for k=1:numel(reqObs)
    if ~isfield(O,reqObs{k})
        error('observations.%s is required.',reqObs{k});
    end
end

reqErr={'sigma_permil','covariance','errorModel'};
for k=1:numel(reqErr)
    if ~isfield(E,reqErr{k})
        error('observationError.%s is required.',reqErr{k});
    end
end

reqDisc={'enabled','observationSD_permil','surrogateCovariance','modelCovariance'};
for k=1:numel(reqDisc)
    if ~isfield(D,reqDisc{k})
        error('modelDiscrepancy.%s is required.',reqDisc{k});
    end
end

obs=O.deltaD_permil(:);
pred=prediction(:);

if isempty(obs) || any(~isfinite(obs))
    error('observations.deltaD_permil must contain finite values.');
end
if numel(pred)~=numel(obs) || any(~isfinite(pred))
    error('prediction must contain one finite value per observation.');
end

hom=O.homologues(:);
if isempty(hom)
    error('observations.homologues must be nonempty.');
end

design=string(O.design);
if strlength(design)==0
    error('observations.design must be nonempty.');
end

switch design
    case "C29_only"
        if ~isequal(hom,29) || numel(obs)~=1
            error('C29_only requires one observation for homologue 29.');
        end

    case "multi_homologue_C27_C29_C31_C33"
        if ~isequal(hom,[27;29;31;33]) || numel(obs)~=4
            error(['multi_homologue_C27_C29_C31_C33 requires four ' ...
                'observations ordered [27 29 31 33].']);
        end

    otherwise
        error('Unsupported observation design.');
end

if ~isempty(E.covariance)
    error('observationError.covariance is not supported.');
end

if string(E.errorModel)~="independent_gaussian"
    error('Only independent_gaussian observation error is supported.');
end

sig=E.sigma_permil(:);
if isempty(sig)
    error('observationError.sigma_permil is required.');
end
if isscalar(sig)
    sig=repmat(sig,numel(obs),1);
end
if numel(sig)~=numel(obs) || any(~isfinite(sig)) || any(sig<=0)
    error(['observationError.sigma_permil must be positive finite scalar ' ...
        'or one value per observation.']);
end

if ~isempty(D.observationSD_permil)
    error(['modelDiscrepancy.observationSD_permil is unsupported by the ' ...
        'accepted N3d2 likelihood semantics and must be empty.']);
end

if ~isempty(D.modelCovariance)
    error(['modelDiscrepancy.modelCovariance is unsupported by the ' ...
        'accepted N3d2 likelihood semantics and must be empty.']);
end

S=[];
if logical(D.enabled)
    S=D.surrogateCovariance;
    if isempty(S)
        S=[];
    else
        if ~isnumeric(S) || ~ismatrix(S) || ...
                size(S,1)~=numel(obs) || size(S,2)~=numel(obs) || ...
                any(~isfinite(S(:)))
            error(['modelDiscrepancy.surrogateCovariance must be a finite ' ...
                'nObservation x nObservation matrix.']);
        end
        tol=1e-12*max(1,norm(S,'fro'));
        if norm(S-S.','fro')>tol
            error('modelDiscrepancy.surrogateCovariance must be symmetric.');
        end
        S=0.5*(S+S.');
        ev=eig(S);
        if min(ev)<-max(1e-12,1e-12*max(abs(ev)))
            error('modelDiscrepancy.surrogateCovariance must be positive semidefinite.');
        end
    end
else
    if ~isempty(D.surrogateCovariance)
        error(['modelDiscrepancy.surrogateCovariance must be empty when ' ...
            'modelDiscrepancy.enabled is false.']);
    end
end

t=O.timeDays(:);
if ~isempty(t)
    if any(~isfinite(t))
        error('observations.timeDays must be finite when provided.');
    end
    if ~(isscalar(t) || numel(t)==numel(obs))
        error('observations.timeDays must be scalar or one value per observation.');
    end
end

if ~isscalar(O.isStandingWax) || ...
        ~(islogical(O.isStandingWax) || isnumeric(O.isStandingWax))
    error('observations.isStandingWax must be a scalar logical value.');
end

L=struct;
L.version="N3_observation_likelihood_schema_production_v1";
L.observedDeltaD_permil=obs;
L.predictedDeltaD_permil=pred;
L.sigma_permil=sig;
L.discrepancyCovariance=S;
L.errorModel="independent_gaussian";
L.observationDesign=design;
L.homologues=hom;
L.timeDays=t;
L.isStandingWax=logical(O.isStandingWax);

L.sourcePolicy=struct;
L.sourcePolicy.observationCovariance="BLOCKED";
L.sourcePolicy.surrogateCovariance="ALLOWED_IF_ENABLED";
L.sourcePolicy.modelCovariance="BLOCKED";
L.sourcePolicy.observationSD_permil="BLOCKED";
L.sourcePolicy.missingObservations="BLOCKED";
L.sourcePolicy.discrepancyComposition= ...
    "Sigma_total = Sigma_observation + Sigma_surrogate";

end
