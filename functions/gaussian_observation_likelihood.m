function [logLikelihood,detail] = gaussian_observation_likelihood( ...
    observedDeltaD_permil,predictedDeltaD_permil,sigmaObs_permil, ...
    observationCovariance_permil2,discrepancyCovariance_permil2)
%GAUSSIAN_OBSERVATION_LIKELIHOOD
% Production observation likelihood kernel.


if nargin < 4 || isempty(observationCovariance_permil2)
    observationCovariance_permil2 = [];
end
if nargin < 5 || isempty(discrepancyCovariance_permil2)
    discrepancyCovariance_permil2 = [];
end

if ~isnumeric(observedDeltaD_permil) || ~isreal(observedDeltaD_permil)
    error('N3c:ObservedType', ...
        'observedDeltaD_permil must be a real numeric vector.');
end
if ~isnumeric(predictedDeltaD_permil) || ~isreal(predictedDeltaD_permil)
    error('N3c:PredictedType', ...
        'predictedDeltaD_permil must be a real numeric vector.');
end
if ~isvector(observedDeltaD_permil) || isempty(observedDeltaD_permil)
    error('N3c:ObservedShape', ...
        'observedDeltaD_permil must be a nonempty vector.');
end
if ~isvector(predictedDeltaD_permil) || isempty(predictedDeltaD_permil)
    error('N3c:PredictedShape', ...
        'predictedDeltaD_permil must be a nonempty vector.');
end

obs = observedDeltaD_permil(:);
pred = predictedDeltaD_permil(:);
n = numel(obs);

if numel(pred) ~= n
    error('N3c:LengthMismatch', ...
        'Observed and predicted vectors must have identical length.');
end
if any(~isfinite(obs))
    error('N3c:NonfiniteObservation', ...
        'All observations must be finite.');
end
if any(~isfinite(pred))
    error('N3c:NonfinitePrediction', ...
        'All predictions must be finite.');
end

if ~isnumeric(sigmaObs_permil) || ~isreal(sigmaObs_permil) || ...
        isempty(sigmaObs_permil) || ~isvector(sigmaObs_permil)
    error('N3c:SigmaType', ...
        'sigmaObs_permil must be a real positive scalar or vector.');
end
sig = sigmaObs_permil(:);
if isscalar(sig)
    sig = repmat(sig,n,1);
elseif numel(sig) ~= n
    error('N3c:SigmaLength', ...
        'Vector sigmaObs_permil must match observation length.');
end
if any(~isfinite(sig)) || any(sig <= 0)
    error('N3c:SigmaSupport', ...
        'sigmaObs_permil must contain only finite positive values.');
end

if ~isempty(observationCovariance_permil2)
    error('N3c:ObservationCovarianceUnsupported', ...
        'Nonempty observation-error covariance.');
end

SigmaObs = diag(sig.^2);
SigmaDisc = zeros(n);

if ~isempty(discrepancyCovariance_permil2)
    C = discrepancyCovariance_permil2;
    if ~isnumeric(C) || ~isreal(C) || ~isequal(size(C),[n n])
        error('N3c:DiscrepancyCovarianceShape', ...
            'Discrepancy covariance must be a real n-by-n numeric matrix.');
    end
    if any(~isfinite(C),'all')
        error('N3c:DiscrepancyCovarianceFinite', ...
            'Discrepancy covariance must be finite.');
    end

    scale = max(1,max(abs(C),[],'all'));
    symTol = 1e-12*scale;
    if max(abs(C-C.'),[],'all') > symTol
        error('N3c:DiscrepancyCovarianceSymmetry', ...
            'Discrepancy covariance must be symmetric.');
    end

    C = 0.5*(C+C.');
    ev = eig(C);
    psdTol = 1e-12*max(1,max(abs(ev)));
    if any(ev < -psdTol)
        error('N3c:DiscrepancyCovariancePSD', ...
            'Discrepancy covariance must be positive semidefinite.');
    end
    SigmaDisc = C;
end

SigmaTotal = SigmaObs + SigmaDisc;
SigmaTotal = 0.5*(SigmaTotal+SigmaTotal.');

[L,p] = chol(SigmaTotal,'lower');
if p ~= 0
    error('N3c:TotalCovariancePD', ...
        'Total likelihood covariance must be positive definite.');
end

residual = obs-pred;
whitenedResidual = L\residual;
quadratic = sum(whitenedResidual.^2);
logDetSigma = 2*sum(log(diag(L)));
normalization = n*log(2*pi) + logDetSigma;

logLikelihood = -0.5*(quadratic + normalization);

detail = struct;
detail.version = "N3_gaussian_observation_likelihood_production_v1";
detail.nObservations = n;
detail.observedDeltaD_permil = obs;
detail.predictedDeltaD_permil = pred;
detail.residual_permil = residual;
detail.sigmaObs_permil = sig;
detail.observationCovariance_permil2 = [];
detail.discrepancyCovariance_permil2 = SigmaDisc;
detail.totalCovariance_permil2 = SigmaTotal;
detail.whitenedResidual = whitenedResidual;
detail.quadraticForm = quadratic;
detail.logDetCovariance = logDetSigma;
detail.gaussianNormalizationTerm = normalization;
detail.logLikelihood = logLikelihood;
detail.usedDiscrepancyCovariance = ~isempty(discrepancyCovariance_permil2);
detail.observationCovarianceBlocked = true;

end
