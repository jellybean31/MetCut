function T = biochemical_latent_mode12_transform(z,options)
%BIOCHEMICAL_LATENT_MODE12_TRANSFORM
%
% Inference transform between biochemical latent Mode 1/2
% coordinates and mechanistic B01-B07 parameters.
%
% INPUT
% -----
% z
%   n x 2 matrix of latent scores
%   finite-flush basis.
%
% options (optional struct)
%   .basisCheckpointPath
%       
%   .priorMean
%       1x2 latent prior mean. Default [0 0].
%   .priorSD
%       1x2 latent prior SD. Default [1 1].
%
% OUTPUT
% ------
% T.B
%   n x 7 reconstructed mechanistic B01-B07 values.
%
% T.feasible
%   n x 1 logical flag. TRUE only when all B01-B07 values lie inside the
%   envelope bounds.
%
% T.logPriorUnnormalized
%   n x 1 unnormalized independent Gaussian latent log-prior, with -Inf
%   assigned to infeasible proposals.


if nargin<2 || isempty(options), options=struct; end
if ~isfield(options,'priorMean'), options.priorMean=[0 0]; end
if ~isfield(options,'priorSD'), options.priorSD=[1 1]; end

z=double(z);

if isvector(z) && numel(z)==2
    z=reshape(z,1,2);
end

if size(z,2)~=2 || any(~isfinite(z),'all')
    error('z must be a finite n x 2 matrix.');
end

priorMean=reshape(double(options.priorMean),1,[]);
priorSD=reshape(double(options.priorSD),1,[]);

if numel(priorMean)~=2 || numel(priorSD)~=2 || any(priorSD<=0)
    error('priorMean and priorSD must each contain two valid values.');
end

if isfield(options,'basisCheckpointPath') && ~isempty(options.basisCheckpointPath)
    cp=char(options.basisCheckpointPath);
else
    cp=resolveCheckpointLocal();
end

if ~isfile(cp)
    error('checkpoint not found: %s',cp);
end

Q=load(cp,'R');
K=Q.R;
F=K.finiteFlush;

bPositions=6:12;

Bmean=reshape(F.Xmean(bPositions),1,[]);
Bsd=reshape(F.Xsd(bPositions),1,[]);
V=F.parameterLoadingsNumeric(:,1:2);

if size(V,1)~=7 || size(V,2)~=2
    error(' latent loading matrix is not 7 x 2.');
end

D=gsa_isotope_parameter_definitions();

Bids=["B01","B02","B03","B04","B05","B06","B07"];
lo=NaN(1,7);
hi=NaN(1,7);

for ib=1:7
    j=find(string({D.id})==Bids(ib),1);
    if isempty(j)
        error('Could not find %s in gsa_isotope_parameter_definitions.',Bids(ib));
    end
    lo(ib)=D(j).low;
    hi(ib)=D(j).high;
end

%% Latent -> mechanistic B transform
Zstd=z*V.';
B=Bmean + Zstd.*Bsd;

feasible=all(B>=lo & B<=hi,2);

%% Unnormalized truncated Gaussian latent prior
u=(z-priorMean)./priorSD;
logPrior=-0.5.*sum(u.^2,2)-sum(log(priorSD));

logPrior(~feasible)=-Inf;

%% Assemble
T=struct;
T.version='2026-09-03_biochemical_latent_mode12_transform_v1';
T.basisCheckpointPath=cp;
T.z=z;
T.B=B;
T.feasible=feasible;
T.logPriorUnnormalized=logPrior;
T.priorMean=priorMean;
T.priorSD=priorSD;
T.Bids=Bids;
T.Bmean=Bmean;
T.Bsd=Bsd;
T.Blow=lo;
T.Bhigh=hi;
T.loadingMatrix=V;
T.standardizedB=Zstd;
T.definition='Mode1/Mode2 latent inference transform with B01-B07 feasibility.';
end

%% =========================================================================
function cp=resolveCheckpointLocal()

p=which('gsa_isotope_parameter_definitions');

if isempty(p)
    error('Cannot resolve project root because gsa_isotope_parameter_definitions is not on MATLAB path.');
end

functionsDir=fileparts(p);
modelRoot=fileparts(functionsDir);

cp=fullfile(modelRoot,'data','reference', ...
    'K3_2E_biochemical_latent_mode_audit_v3.mat');
end
