function [J,detail] = effective7_structural4_linear_observation_objective( ...
    Zvec,M,intracellularNode,timeDay,P,H,yObs,sigmaObs,options)
%EFFECTIVE7_STRUCTURAL4_LINEAR_OBSERVATION_OBJECTIVE
% Gaussian MAP objective for arbitrary linear observations of the four wax
% homolog outputs.
%
% Model:
%   wax = f(Z)                       [4x1]
%   yhat = H * wax                   [m x 1]
%
% Examples:
%   direct four-homolog observations:
%       H = eye(4)
%
%   C29-C27 and C31-C29 spacings:
%       H = [-1 1 0 0;
%             0 -1 1 0]
%
% options.priorWeight default = 1.


if nargin<9 || isempty(options), options=struct; end
if ~isfield(options,'priorWeight') || isempty(options.priorWeight)
    options.priorWeight=1;
end

Z=Zvec(:);
if numel(Z)~=4
    error('Zvec must contain 4 elements.');
end

if size(H,2)~=4
    error('H must have four columns corresponding to C27/C29/C31/C33.');
end

yObs=yObs(:);
if size(H,1)~=numel(yObs)
    error('Number of H rows must match yObs length.');
end

sig=sigmaObs(:);
if isscalar(sig)
    sig=repmat(sig,numel(yObs),1);
end
if numel(sig)~=numel(yObs) || any(sig<=0) || any(~isfinite(sig))
    error('sigmaObs must be positive scalar or match yObs length.');
end

ok=isfinite(yObs);
if ~any(ok)
    error('At least one finite observation is required.');
end

F=effective7_structural4_inverse_forward( ...
    M,intracellularNode,timeDay,P,Z);

wax=F.waxDeltaD_permil(:,1);
yhat=H*wax;

r=(yObs(ok)-yhat(ok))./sig(ok);

dataNLL=0.5*sum(r.^2 + log(2*pi*sig(ok).^2));
priorNLL=0.5*sum(Z.^2);

J=dataNLL + options.priorWeight*priorNLL;

detail=struct;
detail.Z=Z;
detail.waxDeltaD_permil=wax;
detail.H=H;
detail.predictedObservation=yhat;
detail.observedObservation=yObs;
detail.sigmaObservation=sig;
detail.finiteMask=ok;
detail.standardizedResidual=r;
detail.dataNLL=dataNLL;
detail.priorNLL=priorNLL;
detail.priorWeight=options.priorWeight;
detail.objective=J;
detail.forward=F;

end
