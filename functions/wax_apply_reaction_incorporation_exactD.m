function P = wax_apply_reaction_incorporation_exactD(postPools,reactionTraits)
%WAX_APPLY_REACTION_INCORPORATION_EXACTD
% Resolve post-malonate H-source pools into reaction-specific H classes and
% apply reaction-incorporation isotope effects in isotope-ratio space.
%
% INPUT postPools
%   Output of wax_hydrogen_source_pools_postmalonate_exactD.
%
% INPUT reactionTraits
%   Output-compatible structure from wax_reaction_incorporation_traits_null.
%
% Incorporation definition:
%
%   R_incorporated = alpha_reaction .* R_source
%
% where alpha is from:
%   - source-pool alphaFormation
%   - NAD(P)H-water exchange
%   - post-malonate exchange



if ~isstruct(postPools) || ~isfield(postPools,'classes') || ...
        ~isfield(postPools,'sourceClassNames') || ~isfield(postPools,'R_VSMOW')
    error('postPools must be output from wax_hydrogen_source_pools_postmalonate_exactD.');
end
if nargin < 2 || isempty(reactionTraits)
    reactionTraits = wax_reaction_incorporation_traits_null();
end
if ~isstruct(reactionTraits)
    error('reactionTraits must be a structure.');
end

requiredPost = { ...
    'acetylStarterH', ...
    'plastidRetainedMalonylH', ...
    'erRetainedMalonylH', ...
    'plastidNADPHHydrideH', ...
    'plastidNADHHydrideH', ...
    'plastidExchangeableH', ...
    'erNADPHHydrideH', ...
    'erExchangeableH', ...
    'cer3AldehydeHydrideH'};
for i=1:numel(requiredPost)
    if ~isfield(postPools.classes,requiredPost{i})
        error('postPools.classes.%s is required.',requiredPost{i});
    end
end

Rstd = postPools.R_VSMOW;
if isfield(postPools,'timeDays')
    timeDays = postPools.timeDays(:).';
else
    timeDays = 1:numel(postPools.classes.acetylStarterH.ratio);
end
nTime = numel(timeDays);

alpha = struct;
alpha.plastidKARNADPHHydrideH = getAlpha(reactionTraits,{'plastid','KAR','alphaHydride'},nTime);
alpha.plastidENRNADHHydrideH = getAlpha(reactionTraits,{'plastid','ENR','alphaHydride'},nTime);
alpha.plastidENRProtonH = getAlpha(reactionTraits,{'plastid','ENR','alphaProton'},nTime);
alpha.erKCRNADPHHydrideH = getAlpha(reactionTraits,{'er','KCR','alphaHydride'},nTime);
alpha.erECRNADPHHydrideH = getAlpha(reactionTraits,{'er','ECR','alphaHydride'},nTime);
alpha.erECRProtonH = getAlpha(reactionTraits,{'er','ECR','alphaProton'},nTime);
alpha.cer3AldehydeHydrideH = getAlpha(reactionTraits,{'cer3','alphaAldehydeHydride'},nTime);

pre = struct;
pre.acetylStarterH = clonePool(postPools.classes.acetylStarterH);
pre.plastidRetainedMalonylH = clonePool(postPools.classes.plastidRetainedMalonylH);
pre.erRetainedMalonylH = clonePool(postPools.classes.erRetainedMalonylH);
pre.plastidKARNADPHHydrideH = clonePool(postPools.classes.plastidNADPHHydrideH);
pre.plastidENRNADHHydrideH = clonePool(postPools.classes.plastidNADHHydrideH);
pre.plastidENRProtonH = clonePool(postPools.classes.plastidExchangeableH);
pre.erKCRNADPHHydrideH = clonePool(postPools.classes.erNADPHHydrideH);
pre.erECRNADPHHydrideH = clonePool(postPools.classes.erNADPHHydrideH);
pre.erECRProtonH = clonePool(postPools.classes.erExchangeableH);
pre.cer3AldehydeHydrideH = clonePool(postPools.classes.cer3AldehydeHydrideH);

classes = struct;
classes.acetylStarterH = pre.acetylStarterH;
classes.plastidRetainedMalonylH = pre.plastidRetainedMalonylH;
classes.erRetainedMalonylH = pre.erRetainedMalonylH;
classes.plastidKARNADPHHydrideH = applyAlpha(pre.plastidKARNADPHHydrideH,alpha.plastidKARNADPHHydrideH,Rstd,'plastid KAR NADPH hydride');
classes.plastidENRNADHHydrideH = applyAlpha(pre.plastidENRNADHHydrideH,alpha.plastidENRNADHHydrideH,Rstd,'plastid ENR NADH hydride');
classes.plastidENRProtonH = applyAlpha(pre.plastidENRProtonH,alpha.plastidENRProtonH,Rstd,'plastid ENR proton');
classes.erKCRNADPHHydrideH = applyAlpha(pre.erKCRNADPHHydrideH,alpha.erKCRNADPHHydrideH,Rstd,'ER KCR NADPH hydride');
classes.erECRNADPHHydrideH = applyAlpha(pre.erECRNADPHHydrideH,alpha.erECRNADPHHydrideH,Rstd,'ER ECR NADPH hydride');
classes.erECRProtonH = applyAlpha(pre.erECRProtonH,alpha.erECRProtonH,Rstd,'ER ECR proton');
classes.cer3AldehydeHydrideH = applyAlpha(pre.cer3AldehydeHydrideH,alpha.cer3AldehydeHydrideH,Rstd,'CER3 aldehydic hydride');

names = { ...
    'acetylStarterH', ...
    'plastidRetainedMalonylH', ...
    'erRetainedMalonylH', ...
    'plastidKARNADPHHydrideH', ...
    'plastidENRNADHHydrideH', ...
    'plastidENRProtonH', ...
    'erKCRNADPHHydrideH', ...
    'erECRNADPHHydrideH', ...
    'erECRProtonH', ...
    'cer3AldehydeHydrideH'};

P = struct;
P.timeDays = timeDays;
P.R_VSMOW = Rstd;
P.classes = classes;
P.sourceClassNames = names;
P.preIncorporationClasses = pre;
P.alphaByClass = alpha;
P.reactionTraits = reactionTraits;
P.upstreamPostMalonatePools = postPools;
P.sourcePoolFormationAlreadyApplied = true;
P.cofactorWaterExchangeAlreadyApplied = true;
P.postMalonateExchangeAlreadyApplied = true;
P.incorporationFractionationApplied = anyAlphaDiffersFromOne(alpha);
P.isotopeBookkeeping = 'exactDAtomFraction';
P.modelRole = ['reaction-specific incorporated-H isotope states after ', ...
    'source formation/exchange and before whole-wax exact-D combination'];
end

function a = getAlpha(S,path,nTime)
x=S;
for i=1:numel(path)
    if ~isstruct(x) || ~isfield(x,path{i})
        error('Missing reaction-incorporation field: %s',strjoin(path,'.'));
    end
    x=x.(path{i});
end
if isscalar(x)
    a=repmat(x,1,nTime);
else
    a=x(:).';
end
if numel(a)~=nTime || any(~isfinite(a)) || any(a<=0)
    error('%s must be a positive finite scalar or one value per timestep.',strjoin(path,'.'));
end
end

function C = clonePool(C)
% Normalize time-oriented isotope fields without altering upstream metadata.
C.ratio=C.ratio(:).';
C.deltaD=C.deltaD(:).';
C.DAtomFraction=C.DAtomFraction(:).';
end

function C = applyAlpha(C,a,Rstd,label)
r0=C.ratio(:).';
if numel(r0)~=numel(a)
    error('Alpha time dimension does not match source pool for %s.',label);
end
r=a.*r0;
R=Rstd.*r;
x=R./(1+R);
C.preIncorporationRatio=r0;
C.preIncorporationDeltaD=C.deltaD(:).';
C.preIncorporationDAtomFraction=C.DAtomFraction(:).';
C.alphaIncorporation=a;
C.ratio=r;
C.deltaD=1000.*(r-1);
C.DAtomFraction=x;
C.incorporationReaction=label;
end

function tf = anyAlphaDiffersFromOne(A)
fn=fieldnames(A);
tf=false;
for i=1:numel(fn)
    if any(abs(A.(fn{i})-1)>1e-15)
        tf=true;
        return;
    end
end
end
