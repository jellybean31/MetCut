function P = wax_hydrogen_source_pools_postmalonate_exactD(timeDays,intracellular,sourceTraits)
%WAX_HYDROGEN_SOURCE_POOLS_POSTMALONATE_EXACTD
% Exact-D source pools with explicit acetyl-starter, retained-malonyl, and
% CER3 aldehydic-H bookkeeping.
%
% The already validated wax_hydrogen_source_pools_exactD operator is used
% unchanged for plastid NADPH, plastid NADH, ER/cytosolic NADPH, and local
% intracellular-water isotope states. 
%
% Required new fields:
%   sourceTraits.acetylStarterH
%   sourceTraits.retainedMalonylH.plastid
%   sourceTraits.retainedMalonylH.er
%   sourceTraits.cer3AldehydeH
%
% Retained-malonyl configs may additionally specify:
%   .waterExchangeFraction   scalar or one value per timestep, [0,1]
%   .exchangeWaterSource     plastidWater or cytosolWater
%   .alphaWaterExchange      positive scalar
%
% The pre-exchange source can be plastidWater, cytosolWater,
% acetylStarterH, or prescribed.

if ~isstruct(sourceTraits)
    error('sourceTraits must be a structure.');
end

required = {'acetylStarterH','retainedMalonylH','cer3AldehydeH', ...
    'plastidNADPH','plastidNADH','erNADPH'};
for i=1:numel(required)
    if ~isfield(sourceTraits,required{i})
        error('sourceTraits.%s is required.',required{i});
    end
end
if ~isfield(sourceTraits.retainedMalonylH,'plastid') || ...
        ~isfield(sourceTraits.retainedMalonylH,'er')
    error('sourceTraits.retainedMalonylH requires plastid and er fields.');
end

% Build the input needed by the validated pathway-explicit base operator.
B = sourceTraits;
if ~isfield(B,'acetylMalonylH')
    B.acetylMalonylH = sourceTraits.acetylStarterH;
end
if ~isfield(B,'terminalH')
    B.terminalH.mode = 'cytosolWater';
end

Pbase = wax_hydrogen_source_pools_exactD(timeDays,intracellular,B);
Rstd = Pbase.R_VSMOW;
nTime = numel(Pbase.timeDays);
water = Pbase.water;

resolved = struct;
resolved.plastidNADPH = Pbase.resolvedPools.plastidNADPH;
resolved.plastidNADH = Pbase.resolvedPools.plastidNADH;
resolved.erNADPH = Pbase.resolvedPools.erNADPH;

acetyl = resolveEndmemberNew(sourceTraits.acetylStarterH,water,struct,nTime,Rstd,'acetylStarterH');
resolved.acetylStarterH = acetyl;

pMal = resolveEndmemberNew(sourceTraits.retainedMalonylH.plastid,water,resolved,nTime,Rstd,'plastidRetainedMalonylH');
pMal = applyExchangeNew(pMal,sourceTraits.retainedMalonylH.plastid,water,Rstd,'plastidWater');
resolved.plastidRetainedMalonylH = pMal;

erMal = resolveEndmemberNew(sourceTraits.retainedMalonylH.er,water,resolved,nTime,Rstd,'erRetainedMalonylH');
erMal = applyExchangeNew(erMal,sourceTraits.retainedMalonylH.er,water,Rstd,'cytosolWater');
resolved.erRetainedMalonylH = erMal;

cer3 = resolveCER3(sourceTraits.cer3AldehydeH,water,resolved,nTime,Rstd);
resolved.cer3AldehydeH = cer3;

classes = struct;
classes.acetylStarterH = acetyl;
classes.plastidRetainedMalonylH = pMal;
classes.erRetainedMalonylH = erMal;
classes.plastidNADPHHydrideH = resolved.plastidNADPH;
classes.plastidNADHHydrideH = resolved.plastidNADH;
classes.plastidExchangeableH = water.plastid;
classes.erNADPHHydrideH = resolved.erNADPH;
classes.erExchangeableH = water.cytosol;
classes.cer3AldehydeHydrideH = cer3;

P = struct;
P.timeDays = Pbase.timeDays;
P.R_VSMOW = Rstd;
P.water = water;
P.resolvedPools = resolved;
P.classes = classes;
P.sourceClassNames = { ...
    'acetylStarterH', ...
    'plastidRetainedMalonylH', ...
    'erRetainedMalonylH', ...
    'plastidNADPHHydrideH', ...
    'plastidNADHHydrideH', ...
    'plastidExchangeableH', ...
    'erNADPHHydrideH', ...
    'erExchangeableH', ...
    'cer3AldehydeHydrideH'};
P.pathwayModel = Pbase.pathwayModel;
P.baseValidatedSourcePools = Pbase;
P.incorporationFractionationApplied = false;
P.waxMixingApplied = false;
P.postMalonateExchangeApplied = true;
P.isotopeBookkeeping = 'exactDAtomFraction';
P.modelRole = ['explicit post-malonate/CER3 source states before ', ...
    'reaction-specific lipid incorporation isotope effects'];
end

function pool = resolveEndmemberNew(cfg,water,resolved,nTime,Rstd,label)
if ~isstruct(cfg) || ~isfield(cfg,'source')
    error('%s.source is required.',label);
end
src = lower(char(cfg.source));
switch src
    case 'plastidwater'
        r = water.plastid.ratio;
    case 'cytosolwater'
        r = water.cytosol.ratio;
    case 'acetylstarterh'
        if ~isfield(resolved,'acetylStarterH')
            error('%s cannot reference acetylStarterH before it is resolved.',label);
        end
        r = resolved.acetylStarterH.ratio;
    case 'plastidretainedmalonylh'
        if ~isfield(resolved,'plastidRetainedMalonylH')
            error('%s cannot reference plastidRetainedMalonylH before it is resolved.',label);
        end
        r = resolved.plastidRetainedMalonylH.ratio;
    case 'erretainedmalonylh'
        if ~isfield(resolved,'erRetainedMalonylH')
            error('%s cannot reference erRetainedMalonylH before it is resolved.',label);
        end
        r = resolved.erRetainedMalonylH.ratio;
    case 'prescribed'
        r = prescribedRatio(cfg,nTime);
    otherwise
        error('%s has unknown source: %s',label,src);
end

if isfield(cfg,'alphaFormation') && ~isempty(cfg.alphaFormation)
    a = cfg.alphaFormation;
else
    a = 1;
end
if ~isscalar(a) || ~isfinite(a) || a<=0
    error('%s.alphaFormation must be a positive finite scalar.',label);
end
r = a.*r;
pool = makePool(r,Rstd,label);
pool.source = src;
pool.alphaFormation = a;
pool.preExchangeRatio = pool.ratio;
pool.preExchangeDeltaD = pool.deltaD;
pool.preExchangeDAtomFraction = pool.DAtomFraction;
end

function pool = applyExchangeNew(pool,cfg,water,Rstd,defaultWater)
nTime = numel(pool.ratio);
if isfield(cfg,'waterExchangeFraction') && ~isempty(cfg.waterExchangeFraction)
    f = expand(cfg.waterExchangeFraction,nTime,'waterExchangeFraction');
else
    f = zeros(1,nTime);
end
if any(f<0) || any(f>1)
    error('waterExchangeFraction must lie in [0,1].');
end
if isfield(cfg,'exchangeWaterSource') && ~isempty(cfg.exchangeWaterSource)
    src = lower(char(cfg.exchangeWaterSource));
else
    src = lower(defaultWater);
end
switch src
    case 'plastidwater'
        rW = water.plastid.ratio;
    case 'cytosolwater'
        rW = water.cytosol.ratio;
    otherwise
        error('exchangeWaterSource must be plastidWater or cytosolWater.');
end
if isfield(cfg,'alphaWaterExchange') && ~isempty(cfg.alphaWaterExchange)
    a = cfg.alphaWaterExchange;
else
    a = 1;
end
if ~isscalar(a) || ~isfinite(a) || a<=0
    error('alphaWaterExchange must be a positive finite scalar.');
end
xW = ratioToX(a.*rW,Rstd);
x = (1-f).*pool.preExchangeDAtomFraction + f.*xW;
r = xToRatio(x,Rstd);
old = pool;
pool = makePool(r,Rstd,old.description);
pool.source = old.source;
pool.alphaFormation = old.alphaFormation;
pool.preExchangeRatio = old.preExchangeRatio;
pool.preExchangeDeltaD = old.preExchangeDeltaD;
pool.preExchangeDAtomFraction = old.preExchangeDAtomFraction;
pool.waterExchangeFraction = f;
pool.exchangeWaterSource = src;
pool.alphaWaterExchange = a;
end

function pool = resolveCER3(cfg,water,resolved,nTime,Rstd)
if ~isstruct(cfg) || ~isfield(cfg,'mode')
    error('sourceTraits.cer3AldehydeH.mode is required.');
end
mode = lower(char(cfg.mode));
switch mode
    case 'cytosolwater'
        r=water.cytosol.ratio;
    case 'plastidwater'
        r=water.plastid.ratio;
    case 'ernadph'
        r=resolved.erNADPH.ratio;
    case 'plastidnadph'
        r=resolved.plastidNADPH.ratio;
    case 'plastidnadh'
        r=resolved.plastidNADH.ratio;
    case 'acetylstarterh'
        r=resolved.acetylStarterH.ratio;
    case 'plastidretainedmalonylh'
        r=resolved.plastidRetainedMalonylH.ratio;
    case 'erretainedmalonylh'
        r=resolved.erRetainedMalonylH.ratio;
    case 'prescribed'
        r=prescribedRatio(cfg,nTime);
    otherwise
        error('Unknown cer3AldehydeH.mode: %s',mode);
end
pool=makePool(r,Rstd,'CER3 aldehydic H retained into alkane');
pool.mode=mode;
end

function r = prescribedRatio(cfg,nTime)
if isfield(cfg,'ratioNormalized') && ~isempty(cfg.ratioNormalized)
    r=expand(cfg.ratioNormalized,nTime,'ratioNormalized');
elseif isfield(cfg,'deltaD_permil') && ~isempty(cfg.deltaD_permil)
    d=expand(cfg.deltaD_permil,nTime,'deltaD_permil');
    r=1+d./1000;
else
    error('prescribed source requires ratioNormalized or deltaD_permil.');
end
end

function pool=makePool(r,Rstd,description)
r=r(:).'; R=Rstd.*r; x=R./(1+R);
pool=struct('ratio',r,'deltaD',1000.*(r-1),'DAtomFraction',x, ...
    'description',description);
end

function x=ratioToX(r,Rstd)
R=Rstd.*r; x=R./(1+R);
end

function r=xToRatio(x,Rstd)
if any(~isfinite(x)) || any(x<=0) || any(x>=1)
    error('D atom fraction must lie in (0,1).');
end
R=x./(1-x); r=R./Rstd;
end

function v=expand(v,nTime,label)
if isscalar(v), v=repmat(v,1,nTime); else, v=v(:).'; end
if numel(v)~=nTime || any(~isfinite(v))
    error('%s must be scalar or one finite value per timestep.',label);
end
end
