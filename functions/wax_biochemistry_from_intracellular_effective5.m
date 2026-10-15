
function Bio = wax_biochemistry_from_intracellular_effective5( ...
    timeDays,intracellular,cfg,effectiveStates,reactionAlpha)
%WAX_BIOCHEMISTRY_FROM_INTRACELLULAR_EFFECTIVE5
% Five-state reduced biochemical bridge for inverse identifiability tests.
%
% effectiveStates:
%   .acetylStarter_permil
%   .retainedPrecursor_permil
%   .plastidNADPH_permil
%   .plastidNADH_permil
%   .erNADPH_permil


if nargin<3 || isempty(cfg), cfg=struct; end

if ~isfield(cfg,'routeFractions') || isempty(cfg.routeFractions)
    cfg.routeFractions=struct('C16',0,'C18',1);
end

if ~isfield(cfg,'alkaneCarbonNumber') || isempty(cfg.alkaneCarbonNumber)
    cfg.alkaneCarbonNumber=[27 29 31 33];
end

P=wax_integrated_biochemical_prior_definitions();

if nargin<5 || isempty(reactionAlpha)
    reactionAlpha=0.5.*( ...
        P.core.reaction.lowerAlpha(:).' + ...
        P.core.reaction.upperAlpha(:).');
    reactionAlpha(7)=1;
end

reactionAlpha=reactionAlpha(:).';
reactionAlpha(7)=1;

required={ ...
    'acetylStarter_permil', ...
    'retainedPrecursor_permil', ...
    'plastidNADPH_permil', ...
    'plastidNADH_permil', ...
    'erNADPH_permil'};

for k=1:numel(required)
    if ~isfield(effectiveStates,required{k})
        error('effectiveStates.%s is required.',required{k});
    end
end

t=timeDays(:).';
nT=numel(t);

dStarter=expandLocal(effectiveStates.acetylStarter_permil,nT);
dRet=expandLocal(effectiveStates.retainedPrecursor_permil,nT);
dPNP=expandLocal(effectiveStates.plastidNADPH_permil,nT);
dPNH=expandLocal(effectiveStates.plastidNADH_permil,nT);
dERNP=expandLocal(effectiveStates.erNADPH_permil,nT);

T=wax_hydrogen_source_traits_postmalonate_null();

T.acetylStarterH=prescribedLocal(dStarter);
T.acetylMalonylH=T.acetylStarterH;

T.plastidNADPH.pathways.FNR=prescribedLocal(dPNP);
T.plastidNADPH.pathways.OPPP.G6PDH=prescribedLocal(dPNP);
T.plastidNADPH.pathways.OPPP.SixPGDH=prescribedLocal(dPNP);
T.plastidNADPH.pathways.NADPME=prescribedLocal(dPNP);
T.plastidNADPH.waterExchangeFraction=0;
T.plastidNADPH.alphaWaterExchange=1;

T.plastidNADH.pathways.pPDH=prescribedLocal(dPNH);
T.plastidNADH.pathways.GAPDH=prescribedLocal(dPNH);
T.plastidNADH.pathways.MDH=prescribedLocal(dPNH);
T.plastidNADH.pathways.PGDH=prescribedLocal(dPNH);
T.plastidNADH.waterExchangeFraction=0;

T.erNADPH.pathways.OPPP.G6PDH=prescribedLocal(dERNP);
T.erNADPH.pathways.OPPP.SixPGDH=prescribedLocal(dERNP);
T.erNADPH.pathways.ICDH=prescribedLocal(dERNP);
T.erNADPH.pathways.npGAPDH=prescribedLocal(dERNP);
T.erNADPH.pathways.NADPME=prescribedLocal(dERNP);
T.erNADPH.waterExchangeFraction=0;
T.erNADPH.alphaWaterExchange=1;

T.retainedMalonylH.plastid.source='prescribed';
T.retainedMalonylH.plastid.deltaD_permil=dRet;
T.retainedMalonylH.plastid.waterExchangeFraction=0;
T.retainedMalonylH.plastid.alphaWaterExchange=1;

T.retainedMalonylH.er.source='prescribed';
T.retainedMalonylH.er.deltaD_permil=dRet;
T.retainedMalonylH.er.waterExchangeFraction=0;
T.retainedMalonylH.er.alphaWaterExchange=1;

T.cer3AldehydeH.mode=P.cer3Source;

sourcePools=wax_hydrogen_source_pools_postmalonate_exactD(t,intracellular,T);
Rtraits=alphaVectorToTraitsLocal(reactionAlpha);
incorporatedPools=wax_apply_reaction_incorporation_exactD(sourcePools,Rtraits);

S=alkane_biosynthesis_h_stoichiometry_reaction_explicit( ...
    cfg.alkaneCarbonNumber,cfg.routeFractions);

W=wax_combine_hydrogen_classes_exactD(S,incorporatedPools);

Bio=struct;
Bio.version='2026-08-28_effective5_inverse_bridge_v1';
Bio.timeDays=t;
Bio.parameterMode='effective5InverseDiagnostic';
Bio.effectiveStates=effectiveStates;
Bio.sourceTraits=T;
Bio.sourcePools=sourcePools;
Bio.reactionTraits=Rtraits;
Bio.incorporatedPools=incorporatedPools;
Bio.stoichiometry=S;
Bio.alkaneCarbonNumber=cfg.alkaneCarbonNumber(:).';
Bio.wax=W;

end

function E=prescribedLocal(d)
E=struct('source','prescribed','deltaD_permil',d, ...
    'ratioNormalized',[],'alphaFormation',1);
end

function x=expandLocal(x,n)
if isscalar(x)
    x=repmat(x,1,n);
else
    x=x(:).';
    if numel(x)~=n
        error('Effective-state series must be scalar or match timeDays.');
    end
end
end

function T=alphaVectorToTraitsLocal(X)
T=wax_reaction_incorporation_traits_null();
T.plastid.KAR.alphaHydride=X(1);
T.plastid.ENR.alphaHydride=X(2);
T.plastid.ENR.alphaProton=X(3);
T.er.KCR.alphaHydride=X(4);
T.er.ECR.alphaHydride=X(5);
T.er.ECR.alphaProton=X(6);
T.cer3.alphaAldehydeHydride=X(7);
end
