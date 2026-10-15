function P = wax_hydrogen_source_pools_exactD(timeDays,intracellular,sourceTraits)
%WAX_HYDROGEN_SOURCE_POOLS_EXACTD
%
% Pathway-explicit exact-D isotope-state operator for H source pools used
% by the plant n-alkane biosynthesis stoichiometry model.
%
% This version replaces the former generic arbitrary "components" model
% with named metabolic pathway pools:
%
%   inherited organic precursor H
%       -> acetyl/malonyl H pool
%
%   plastid NADPH
%       -> FNR
%       -> OPPP
%       -> NADP-ME
%
%   plastid NADH
%       -> pPDH
%       -> GAPDH
%       -> MDH
%       -> PGDH
%
%   ER/cytosolic NADPH
%       -> OPPP
%       -> ICDH
%       -> npGAPDH
%       -> NADP-ME
%
% -------------------------------------------------------------------------
% REQUIRED INTRACELLULAR STATES
%
% intracellular.ratio.plastid or intracellular.deltaD.plastid
% intracellular.ratio.cytosol or intracellular.deltaD.cytosol
%
% -------------------------------------------------------------------------
% REQUIRED sourceTraits STRUCTURE
%
% sourceTraits.acetylMalonylH
% sourceTraits.plastidNADPH
% sourceTraits.plastidNADH
% sourceTraits.erNADPH
% sourceTraits.terminalH
%
% See wax_hydrogen_source_traits_pathway_null.m for a complete template.
%
% -------------------------------------------------------------------------
% ISOTOPE END-MEMBER SPECIFICATION
%
% Any pathway H end-member can be defined by one of:
%
%   .source = 'plastidWater'
%   .source = 'cytosolWater'
%   .source = 'acetylMalonylH'
%   .source = 'prescribed'
%
% For prescribed:
%
%   .deltaD_permil
%
% or
%
%   .ratioNormalized
%
% Scalar or one value per timestep is accepted.
%
% Each end-member may additionally specify:
%
%   .alphaFormation
%
% with:
%
%   R_formed = alphaFormation * R_source
%
% alphaFormation belongs to SOURCE-POOL FORMATION. 
%
% -------------------------------------------------------------------------
% PATHWAY FRACTIONS
%
% Fractions can be scalar or one value per timestep.
% They must be nonnegative and sum to 1 at every timestep.
%
% sourceTraits.plastidNADPH.pathwayFractions:
%   .FNR
%   .OPPP
%   .NADPME
%
% sourceTraits.plastidNADH.pathwayFractions:
%   .pPDH
%   .GAPDH
%   .MDH
%   .PGDH
%
% sourceTraits.erNADPH.pathwayFractions:
%   .OPPP
%   .ICDH
%   .npGAPDH
%   .NADPME
%
% -------------------------------------------------------------------------
% OPPP STOICHIOMETRY
%
% Each complete OPPP oxidative pass produces two NADPH:
%
%   1 from G6PDH
%   1 from 6PGDH
%
% Therefore the OPPP active-hydride end-member is internally:
%
%   x_OPPP = 0.5*x_G6PDH + 0.5*x_6PGDH
%
% in exact D atom-fraction space.
%
% This fixed 50:50 ratio is reaction stoichiometry, not a tunable pathway
% fraction.
%
% -------------------------------------------------------------------------
% EFFECTIVE ACTIVE-H WATER EXCHANGE
%
% NADPH/NADH pools can optionally undergo effective exchange with local
% water after formation:
%
%   x_final =
%       (1-f_exchange)*x_formed
%       + f_exchange*x_exchangeTarget
%
% Fields:
%
%   .waterExchangeFraction
%   .exchangeWaterSource
%   .alphaWaterExchange
%
% Default exchange fraction is 0 if omitted.
%
% -------------------------------------------------------------------------
% TERMINAL H MODES
%
% sourceTraits.terminalH.mode:
%
%   'cytosolWater'
%   'plastidWater'
%   'erNADPH'
%   'plastidNADPH'
%   'plastidNADH'
%   'acetylMalonylH'
%   'prescribed'
%
% For prescribed supply deltaD_permil or ratioNormalized.
%
% -------------------------------------------------------------------------
% OUTPUT SOURCE CLASSES
%
% P.classes.inheritedOrganicPrecursorH
% P.classes.plastidNADPHHydrideH
% P.classes.plastidNADHHydrideH
% P.classes.plastidExchangeableH
% P.classes.erNADPHHydrideH
% P.classes.erExchangeableH
% P.classes.terminalAlkaneFormingH
%
% -------------------------------------------------------------------------

%% ========================================================================
% Time / standard
% ========================================================================

t = timeDays(:);

if isempty(t) || any(~isfinite(t))
    error('timeDays must contain finite values.');
end

nTime = numel(t);

if isfield(sourceTraits,'R_VSMOW')
    Rstd = sourceTraits.R_VSMOW;
else
    Rstd = 155.76e-6;
end

if ~isscalar(Rstd) || ~isfinite(Rstd) || Rstd <= 0
    error('sourceTraits.R_VSMOW must be a positive finite scalar.');
end

%% ========================================================================
% Reject superseded generic-components configuration
% ========================================================================

legacyNames = {'organicPrecursor'};

for i = 1:numel(legacyNames)
    if isfield(sourceTraits,legacyNames{i})
        error([ ...
            'Legacy generic sourceTraits.%s configuration detected. ', ...
            'Use the pathway-explicit sourceTraits.acetylMalonylH, ', ...
            'plastidNADPH, plastidNADH, erNADPH, and terminalH fields.'], ...
            legacyNames{i});
    end
end

%% ========================================================================
% Intracellular water isotope states
% ========================================================================

rPlastid = selectIntracellularRatio( ...
    intracellular,'plastid',nTime);

rCytosol = selectIntracellularRatio( ...
    intracellular,'cytosol',nTime);

water = struct;

water.plastid = makePoolFromRatio( ...
    rPlastid,Rstd,'intracellular plastid water');

water.cytosol = makePoolFromRatio( ...
    rCytosol,Rstd,'intracellular cytosolic water');

%% ========================================================================
% Required source definitions
% ========================================================================

required = { ...
    'acetylMalonylH', ...
    'plastidNADPH', ...
    'plastidNADH', ...
    'erNADPH', ...
    'terminalH'};

for i = 1:numel(required)

    if ~isfield(sourceTraits,required{i})
        error('sourceTraits.%s is required.',required{i});
    end

end

%% ========================================================================
% 1. Acetyl/malonyl inherited-H pool
% ========================================================================

acetylMalonyl = resolveEndmember( ...
    sourceTraits.acetylMalonylH, ...
    water, ...
    struct, ...
    nTime, ...
    Rstd, ...
    'acetylMalonylH');

%% ========================================================================
% 2. Plastid NADPH pathway mixture
% ========================================================================

resolvedBase = struct;
resolvedBase.acetylMalonylH = acetylMalonyl;

plastidNADPH = resolvePlastidNADPH( ...
    sourceTraits.plastidNADPH, ...
    water, ...
    resolvedBase, ...
    nTime, ...
    Rstd);

%% ========================================================================
% 3. Plastid NADH pathway mixture
% ========================================================================

plastidNADH = resolvePlastidNADH( ...
    sourceTraits.plastidNADH, ...
    water, ...
    resolvedBase, ...
    nTime, ...
    Rstd);

%% ========================================================================
% 4. ER/cytosolic NADPH pathway mixture
% ========================================================================

erNADPH = resolveERNADPH( ...
    sourceTraits.erNADPH, ...
    water, ...
    resolvedBase, ...
    nTime, ...
    Rstd);

%% ========================================================================
% 5. Terminal alkane-forming H
% ========================================================================

resolvedAll = resolvedBase;
resolvedAll.plastidNADPH = plastidNADPH;
resolvedAll.plastidNADH = plastidNADH;
resolvedAll.erNADPH = erNADPH;

terminalH = resolveTerminalH( ...
    sourceTraits.terminalH, ...
    water, ...
    resolvedAll, ...
    nTime, ...
    Rstd);

%% ========================================================================
% Map onto stoichiometry source classes
% ========================================================================

classes = struct;

classes.inheritedOrganicPrecursorH = acetylMalonyl;
classes.plastidNADPHHydrideH = plastidNADPH;
classes.plastidNADHHydrideH = plastidNADH;

classes.plastidExchangeableH = water.plastid;

classes.erNADPHHydrideH = erNADPH;
classes.erExchangeableH = water.cytosol;

classes.terminalAlkaneFormingH = terminalH;

%% ========================================================================
% Output
% ========================================================================

P = struct;

P.timeDays = t.';
P.R_VSMOW = Rstd;

P.water = water;

P.resolvedPools = struct;
P.resolvedPools.acetylMalonylH = acetylMalonyl;
P.resolvedPools.plastidNADPH = plastidNADPH;
P.resolvedPools.plastidNADH = plastidNADH;
P.resolvedPools.erNADPH = erNADPH;
P.resolvedPools.terminalH = terminalH;

P.classes = classes;

P.sourceClassNames = { ...
    'inheritedOrganicPrecursorH', ...
    'plastidNADPHHydrideH', ...
    'plastidNADHHydrideH', ...
    'plastidExchangeableH', ...
    'erNADPHHydrideH', ...
    'erExchangeableH', ...
    'terminalAlkaneFormingH'};

P.pathwayModel = struct;

P.pathwayModel.plastidNADPH = ...
    {'FNR','OPPP','NADPME'};

P.pathwayModel.plastidNADH = ...
    {'pPDH','GAPDH','MDH','PGDH'};

P.pathwayModel.erNADPH = ...
    {'OPPP','ICDH','npGAPDH','NADPME'};

P.incorporationFractionationApplied = false;
P.waxMixingApplied = false;

P.isotopeBookkeeping = 'exactDAtomFraction';

P.modelRole = ...
    ['pathway-explicit isotope states for stoichiometric H-source classes ', ...
     'before reaction-specific lipid incorporation effects'];

end

%% =========================================================================
function pool = resolvePlastidNADPH(cfg,water,resolved,nTime,Rstd)

pathNames = {'FNR','OPPP','NADPME'};

fractions = getPathwayFractions( ...
    cfg,pathNames,nTime,'plastidNADPH');

if ~isfield(cfg,'pathways') || ~isstruct(cfg.pathways)
    error('sourceTraits.plastidNADPH.pathways must be a structure.');
end

% FNR
requirePathway(cfg.pathways,'FNR','plastidNADPH');

FNR = resolveEndmember( ...
    cfg.pathways.FNR, ...
    water,resolved,nTime,Rstd, ...
    'plastidNADPH.FNR');

% OPPP = fixed 50:50 active-H mixture from G6PDH and 6PGDH
requirePathway(cfg.pathways,'OPPP','plastidNADPH');

OPPP = resolveOPPP( ...
    cfg.pathways.OPPP, ...
    water,resolved,nTime,Rstd, ...
    'plastidNADPH.OPPP');

% NADP-ME
requirePathway(cfg.pathways,'NADPME','plastidNADPH');

NADPME = resolveEndmember( ...
    cfg.pathways.NADPME, ...
    water,resolved,nTime,Rstd, ...
    'plastidNADPH.NADPME');

pathPools = {FNR,OPPP,NADPME};

pool = mixPathwayPools( ...
    pathNames,pathPools,fractions,Rstd, ...
    'plastid NADPH active hydride');

pool = applyWaterExchange( ...
    pool,cfg,water,Rstd,'plastidWater');

end

%% =========================================================================
function pool = resolvePlastidNADH(cfg,water,resolved,nTime,Rstd)

pathNames = {'pPDH','GAPDH','MDH','PGDH'};

fractions = getPathwayFractions( ...
    cfg,pathNames,nTime,'plastidNADH');

if ~isfield(cfg,'pathways') || ~isstruct(cfg.pathways)
    error('sourceTraits.plastidNADH.pathways must be a structure.');
end

pathPools = cell(1,numel(pathNames));

for i = 1:numel(pathNames)

    nm = pathNames{i};

    requirePathway(cfg.pathways,nm,'plastidNADH');

    pathPools{i} = resolveEndmember( ...
        cfg.pathways.(nm), ...
        water,resolved,nTime,Rstd, ...
        ['plastidNADH.' nm]);

end

pool = mixPathwayPools( ...
    pathNames,pathPools,fractions,Rstd, ...
    'plastid NADH active hydride');

pool = applyWaterExchange( ...
    pool,cfg,water,Rstd,'plastidWater');

end

%% =========================================================================
function pool = resolveERNADPH(cfg,water,resolved,nTime,Rstd)

pathNames = {'OPPP','ICDH','npGAPDH','NADPME'};

fractions = getPathwayFractions( ...
    cfg,pathNames,nTime,'erNADPH');

if ~isfield(cfg,'pathways') || ~isstruct(cfg.pathways)
    error('sourceTraits.erNADPH.pathways must be a structure.');
end

pathPools = cell(1,numel(pathNames));

for i = 1:numel(pathNames)

    nm = pathNames{i};

    requirePathway(cfg.pathways,nm,'erNADPH');

    if strcmp(nm,'OPPP')

        pathPools{i} = resolveOPPP( ...
            cfg.pathways.OPPP, ...
            water,resolved,nTime,Rstd, ...
            'erNADPH.OPPP');

    else

        pathPools{i} = resolveEndmember( ...
            cfg.pathways.(nm), ...
            water,resolved,nTime,Rstd, ...
            ['erNADPH.' nm]);

    end

end

pool = mixPathwayPools( ...
    pathNames,pathPools,fractions,Rstd, ...
    'ER/cytosolic NADPH active hydride');

pool = applyWaterExchange( ...
    pool,cfg,water,Rstd,'cytosolWater');

end

%% =========================================================================
function pool = resolveOPPP(cfg,water,resolved,nTime,Rstd,label)

if ~isstruct(cfg)
    error('%s configuration must be a structure.',label);
end

if ~isfield(cfg,'G6PDH') || ~isfield(cfg,'SixPGDH')
    error('%s requires G6PDH and SixPGDH end-members.',label);
end

G6 = resolveEndmember( ...
    cfg.G6PDH, ...
    water,resolved,nTime,Rstd, ...
    [label '.G6PDH']);

P6 = resolveEndmember( ...
    cfg.SixPGDH, ...
    water,resolved,nTime,Rstd, ...
    [label '.SixPGDH']);

x = ...
    0.5.*G6.DAtomFraction ...
    + 0.5.*P6.DAtomFraction;

r = atomFractionToRatio(x,Rstd);

pool = makePoolFromRatio( ...
    r,Rstd,[label ' 50:50 G6PDH/6PGDH hydride']);

pool.subpathways = struct;
pool.subpathways.G6PDH = G6;
pool.subpathways.SixPGDH = P6;

pool.subpathwayFractions = struct;
pool.subpathwayFractions.G6PDH = 0.5;
pool.subpathwayFractions.SixPGDH = 0.5;

end

%% =========================================================================
function fractions = getPathwayFractions(cfg,pathNames,nTime,label)

if ~isfield(cfg,'pathwayFractions') || ~isstruct(cfg.pathwayFractions)
    error('sourceTraits.%s.pathwayFractions must be a structure.',label);
end

fractions = zeros(numel(pathNames),nTime);

for i = 1:numel(pathNames)

    nm = pathNames{i};

    if ~isfield(cfg.pathwayFractions,nm)
        error('sourceTraits.%s.pathwayFractions.%s is required.',label,nm);
    end

    fractions(i,:) = expandSeries( ...
        cfg.pathwayFractions.(nm), ...
        nTime, ...
        [label '.pathwayFractions.' nm]);

    if any(fractions(i,:) < 0)
        error('%s pathway fractions must be nonnegative.',label);
    end

end

sumF = sum(fractions,1);

if max(abs(sumF-1)) > 1e-12
    error('%s pathway fractions must sum to 1 at every timestep.',label);
end

end

%% =========================================================================
function pool = mixPathwayPools(pathNames,pathPools,fractions,Rstd,description)

nPath = numel(pathNames);
nTime = numel(pathPools{1}.DAtomFraction);

xPath = zeros(nPath,nTime);

for i = 1:nPath
    xPath(i,:) = pathPools{i}.DAtomFraction;
end

xMix = sum(fractions.*xPath,1);

r = atomFractionToRatio(xMix,Rstd);

pool = makePoolFromRatio(r,Rstd,description);

pool.preExchangeRatio = pool.ratio;
pool.preExchangeDeltaD = pool.deltaD;
pool.preExchangeDAtomFraction = pool.DAtomFraction;

pool.pathwayNames = pathNames;
pool.pathwayFractions = fractions;
pool.pathwayPools = struct;

for i = 1:nPath
    pool.pathwayPools.(pathNames{i}) = pathPools{i};
end

end

%% =========================================================================
function pool = applyWaterExchange(pool,cfg,water,Rstd,defaultSource)

if isfield(cfg,'waterExchangeFraction') && ...
        ~isempty(cfg.waterExchangeFraction)

    fExchange = expandSeries( ...
        cfg.waterExchangeFraction, ...
        numel(pool.DAtomFraction), ...
        'waterExchangeFraction');

else

    fExchange = zeros(size(pool.DAtomFraction));

end

if any(fExchange < 0) || any(fExchange > 1)
    error('waterExchangeFraction must lie in [0,1].');
end

if isfield(cfg,'exchangeWaterSource') && ...
        ~isempty(cfg.exchangeWaterSource)

    waterSource = lower(char(cfg.exchangeWaterSource));

else

    waterSource = lower(defaultSource);

end

switch waterSource

    case 'plastidwater'
        rWater = water.plastid.ratio;

    case 'cytosolwater'
        rWater = water.cytosol.ratio;

    otherwise
        error('exchangeWaterSource must be plastidWater or cytosolWater.');
end

if isfield(cfg,'alphaWaterExchange') && ...
        ~isempty(cfg.alphaWaterExchange)

    alpha = cfg.alphaWaterExchange;

else

    alpha = 1;

end

if ~isscalar(alpha) || ~isfinite(alpha) || alpha <= 0
    error('alphaWaterExchange must be a positive finite scalar.');
end

rTarget = alpha.*rWater;

xTarget = ratioToAtomFraction(rTarget,Rstd);

xFinal = ...
    (1-fExchange).*pool.preExchangeDAtomFraction ...
    + fExchange.*xTarget;

rFinal = atomFractionToRatio(xFinal,Rstd);

finalPool = makePoolFromRatio( ...
    rFinal,Rstd,pool.description);

% Preserve formation/pathway diagnostics.
fields = fieldnames(pool);

for i = 1:numel(fields)

    f = fields{i};

    if ~isfield(finalPool,f)
        finalPool.(f) = pool.(f);
    end

end

finalPool.preExchangeRatio = pool.preExchangeRatio;
finalPool.preExchangeDeltaD = pool.preExchangeDeltaD;
finalPool.preExchangeDAtomFraction = pool.preExchangeDAtomFraction;

finalPool.waterExchangeFraction = fExchange;
finalPool.exchangeWaterSource = waterSource;
finalPool.alphaWaterExchange = alpha;

pool = finalPool;

end

%% =========================================================================
function pool = resolveTerminalH(cfg,water,resolved,nTime,Rstd)

if ~isstruct(cfg) || ~isfield(cfg,'mode')
    error('sourceTraits.terminalH.mode is required.');
end

mode = lower(char(cfg.mode));

switch mode

    case 'cytosolwater'
        r = water.cytosol.ratio;
        description = 'terminal H from cytosolic water';

    case 'plastidwater'
        r = water.plastid.ratio;
        description = 'terminal H from plastid water';

    case 'ernadph'
        r = resolved.erNADPH.ratio;
        description = 'terminal H from ER/cytosolic NADPH';

    case 'plastidnadph'
        r = resolved.plastidNADPH.ratio;
        description = 'terminal H from plastid NADPH';

    case 'plastidnadh'
        r = resolved.plastidNADH.ratio;
        description = 'terminal H from plastid NADH';

    case 'acetylmalonylh'
        r = resolved.acetylMalonylH.ratio;
        description = 'terminal H from acetyl/malonyl H pool';

    case 'prescribed'

        r = resolvePrescribedRatio(cfg,nTime);
        description = 'prescribed terminal H';

    otherwise

        error('Unknown terminalH.mode: %s',mode);

end

pool = makePoolFromRatio(r,Rstd,description);
pool.mode = mode;

end

%% =========================================================================
function pool = resolveEndmember(cfg,water,resolved,nTime,Rstd,label)

if ~isstruct(cfg) || ~isfield(cfg,'source')
    error('%s.source is required.',label);
end

sourceName = lower(char(cfg.source));

switch sourceName

    case 'plastidwater'
        rSource = water.plastid.ratio;

    case 'cytosolwater'
        rSource = water.cytosol.ratio;

    case 'acetylmalonylh'

        if ~isfield(resolved,'acetylMalonylH')
            error('%s cannot reference acetylMalonylH before it is resolved.',label);
        end

        rSource = resolved.acetylMalonylH.ratio;

    case 'prescribed'
        rSource = resolvePrescribedRatio(cfg,nTime);

    otherwise
        error('%s has unknown source: %s',label,sourceName);

end

if isfield(cfg,'alphaFormation') && ~isempty(cfg.alphaFormation)
    alpha = cfg.alphaFormation;
else
    alpha = 1;
end

if ~isscalar(alpha) || ~isfinite(alpha) || alpha <= 0
    error('%s.alphaFormation must be a positive finite scalar.',label);
end

rFormed = alpha.*rSource;

if any(~isfinite(rFormed)) || any(rFormed <= 0)
    error('%s formation produced a nonphysical isotope ratio.',label);
end

pool = makePoolFromRatio( ...
    rFormed,Rstd,label);

pool.source = sourceName;
pool.alphaFormation = alpha;
pool.sourceRatio = rSource;
pool.sourceDeltaD = 1000.*(rSource-1);

end

%% =========================================================================
function r = resolvePrescribedRatio(cfg,nTime)

if isfield(cfg,'ratioNormalized') && ~isempty(cfg.ratioNormalized)

    r = expandSeries( ...
        cfg.ratioNormalized,nTime,'ratioNormalized');

elseif isfield(cfg,'deltaD_permil') && ~isempty(cfg.deltaD_permil)

    deltaD = expandSeries( ...
        cfg.deltaD_permil,nTime,'deltaD_permil');

    r = 1 + deltaD./1000;

else

    error('prescribed source requires ratioNormalized or deltaD_permil.');

end

end

%% =========================================================================
function requirePathway(S,name,parent)

if ~isfield(S,name)
    error('sourceTraits.%s.pathways.%s is required.',parent,name);
end

end

%% =========================================================================
function pool = makePoolFromRatio(r,Rstd,description)

r = r(:).';

R = Rstd.*r;
x = R./(1+R);

pool = struct;

pool.ratio = r;
pool.deltaD = 1000.*(r-1);
pool.DAtomFraction = x;
pool.description = description;

end

%% =========================================================================
function x = ratioToAtomFraction(r,Rstd)

R = Rstd.*r;
x = R./(1+R);

end

%% =========================================================================
function r = atomFractionToRatio(x,Rstd)

if any(~isfinite(x)) || any(x <= 0) || any(x >= 1)
    error('D atom fraction must lie in (0,1).');
end

R = x./(1-x);
r = R./Rstd;

end

%% =========================================================================
function r = selectIntracellularRatio(I,name,nTime)

if isfield(I,'ratio') && ...
        isstruct(I.ratio) && ...
        isfield(I.ratio,name)

    r = I.ratio.(name)(:);

elseif isfield(I,'deltaD') && ...
        isstruct(I.deltaD) && ...
        isfield(I.deltaD,name)

    r = 1 + I.deltaD.(name)(:)./1000;

else

    error('Intracellular isotope state for %s was not found.',name);

end

if numel(r) ~= nTime
    error('Intracellular %s isotope state must match timeDays.',name);
end

r = r(:).';

end

%% =========================================================================
function v = expandSeries(v,nTime,label)

if isscalar(v)
    v = repmat(v,1,nTime);
else
    v = v(:).';
end

if numel(v) ~= nTime || any(~isfinite(v))
    error('%s must be scalar or contain one finite value per timestep.',label);
end

end
