
function T = full_chain_wax_only_truth(Rpilot,mode,seed)
%FULL_CHAIN_WAX_ONLY_TRUTH Build noisy downstream wax observations.

if nargin<2 || isempty(mode), mode="timeResolved"; end
if nargin<3 || isempty(seed), seed=20260828; end

M=Rpilot.truth.model;
Y=M.biochemistry.wax.deltaDMatrix;
C=M.biochemistry.alkaneCarbonNumber(:);

if size(Y,1)~=numel(C) && size(Y,2)==numel(C)
    Y=Y.';
end

target=[27;29;31;33];
idx=zeros(4,1);
for k=1:4
    idx(k)=find(C==target(k),1);
end
Y=Y(idx,:);

t=M.biochemistry.wax.timeDays(:).';

switch string(mode)
    case "final"
        sampleIndex=numel(t);
    case "timeResolved"
        dt=median(diff(t));
        stride=max(1,round(0.5/dt));
        sampleIndex=unique([1:stride:numel(t),numel(t)]);
    otherwise
        error('Unknown mode.');
end

rng(seed,'twister');
sigma=2;
Ytrue=Y(:,sampleIndex);
Yobs=Ytrue+sigma*randn(size(Ytrue));

T=struct;
T.mode=string(mode);
T.seed=seed;
T.sigma_permil=sigma;
T.timeDays=t;
T.sampleIndex=sampleIndex;
T.sampleTimeDays=t(sampleIndex);
T.true_permil=Ytrue;
T.observed_permil=Yobs;
T.homologs=target;
end
