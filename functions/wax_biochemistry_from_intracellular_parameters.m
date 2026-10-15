function Bio = wax_biochemistry_from_intracellular_parameters( ...
    timeDays,intracellular,cfg,sourceVals,reactionAlpha)
%WAX_BIOCHEMISTRY_FROM_INTRACELLULAR_PARAMETERS
% Parameterized exact-D biochemical bridge for GSA/inversion development.
%
% sourceVals is the eight-element source parameter vector used by
% wax_integrated_biochemical_prior_definitions:
%   1 acetyl/malonyl inherited H offset
%   2 plastid NADPH pre-exchange offset
%   3 plastid NADH pre-exchange offset
%   4 ER/cytosolic NADPH pre-exchange offset
%   5 plastid NADPH water-exchange fraction
%   6 ER NADPH water-exchange fraction
%   7 post-malonate plastid exchange fraction
%   8 post-malonate ER exchange fraction
%
% reactionAlpha is the seven-element reaction alpha vector:
%   pKAR, pENR hydride, pENR proton, ER KCR, ER ECR hydride,
%   ER ECR proton, CER3.


if nargin<3 || isempty(cfg), cfg=struct; end
if ~isfield(cfg,'routeFractions') || isempty(cfg.routeFractions)
    cfg.routeFractions=struct('C16',0,'C18',1);
end
if ~isfield(cfg,'alkaneCarbonNumber') || isempty(cfg.alkaneCarbonNumber)
    cfg.alkaneCarbonNumber=[27 29 31 33];
end

if ~isfield(cfg,'architectureVersion') || isempty(cfg.architectureVersion)
    cfg.architectureVersion='v1';
end
if ~isfield(cfg,'erStageStates') || isempty(cfg.erStageStates)
    cfg.erStageStates=wax_er_stage_state_traits();
end

P=wax_integrated_biochemical_prior_definitions();

if nargin<4 || isempty(sourceVals)
    sourceVals=0.5.*(P.core.sourceLower+P.core.sourceUpper);
end
if nargin<5 || isempty(reactionAlpha)
    reactionAlpha=0.5.*( ...
        P.core.reaction.lowerAlpha(:).' + ...
        P.core.reaction.upperAlpha(:).');
    reactionAlpha(7)=1;
end

sourceVals=sourceVals(:).';
reactionAlpha=reactionAlpha(:).';

if numel(sourceVals)~=8
    error('sourceVals must contain eight parameters.');
end
if numel(reactionAlpha)~=7
    error('reactionAlpha must contain seven values.');
end
if any(~isfinite(sourceVals)) || any(~isfinite(reactionAlpha))
    error('Biochemical parameter vectors must be finite.');
end
if any(sourceVals(5:8)<0 | sourceVals(5:8)>1)
    error('Biochemical exchange fractions must lie in [0,1].');
end
if any(reactionAlpha<=0)
    error('Reaction alpha values must be positive.');
end

% CER3 remains fixed at unity in the current accepted architecture.
reactionAlpha(7)=1;

t=timeDays(:).';
nTime=numel(t);

if ~isfield(intracellular,'deltaD') || ...
        ~isfield(intracellular.deltaD,'cytosol') || ...
        ~isfield(intracellular.deltaD,'plastid')
    error('intracellular.deltaD.cytosol and .plastid are required.');
end

dCy=intracellular.deltaD.cytosol(:).';
dPl=intracellular.deltaD.plastid(:).';
if numel(dCy)~=nTime || numel(dPl)~=nTime
    error('Intracellular isotope time series must match timeDays.');
end

refWater=P.sourceReferenceWater_permil;

dAc=wax_source_reference_offset_to_deltaD(sourceVals(1),dPl,refWater);
dPNP=wax_source_reference_offset_to_deltaD(sourceVals(2),dPl,refWater);
dPNH=wax_source_reference_offset_to_deltaD(sourceVals(3),dPl,refWater);
dERNP=wax_source_reference_offset_to_deltaD(sourceVals(4),dCy,refWater);

T=wax_hydrogen_source_traits_postmalonate_null();
T.acetylStarterH=prescribedEndmemberLocal(dAc);
T.acetylMalonylH=T.acetylStarterH;

T.plastidNADPH.pathways.FNR=prescribedEndmemberLocal(dPNP);
T.plastidNADPH.pathways.OPPP.G6PDH=prescribedEndmemberLocal(dPNP);
T.plastidNADPH.pathways.OPPP.SixPGDH=prescribedEndmemberLocal(dPNP);
T.plastidNADPH.pathways.NADPME=prescribedEndmemberLocal(dPNP);

T.plastidNADH.pathways.pPDH=prescribedEndmemberLocal(dPNH);
T.plastidNADH.pathways.GAPDH=prescribedEndmemberLocal(dPNH);
T.plastidNADH.pathways.MDH=prescribedEndmemberLocal(dPNH);
T.plastidNADH.pathways.PGDH=prescribedEndmemberLocal(dPNH);

T.erNADPH.pathways.OPPP.G6PDH=prescribedEndmemberLocal(dERNP);
T.erNADPH.pathways.OPPP.SixPGDH=prescribedEndmemberLocal(dERNP);
T.erNADPH.pathways.ICDH=prescribedEndmemberLocal(dERNP);
T.erNADPH.pathways.npGAPDH=prescribedEndmemberLocal(dERNP);
T.erNADPH.pathways.NADPME=prescribedEndmemberLocal(dERNP);

T.plastidNADPH.waterExchangeFraction=sourceVals(5);
T.plastidNADPH.alphaWaterExchange=1;
T.plastidNADPH.exchangeWaterSource='plastidWater';

T.plastidNADH.waterExchangeFraction=P.plastidNADHWaterExchange.fixedFraction;

T.erNADPH.waterExchangeFraction=sourceVals(6);
T.erNADPH.alphaWaterExchange=1;
T.erNADPH.exchangeWaterSource='cytosolWater';

T.retainedMalonylH.plastid.source='acetylStarterH';
T.retainedMalonylH.plastid.waterExchangeFraction=sourceVals(7);
T.retainedMalonylH.plastid.exchangeWaterSource='plastidWater';
T.retainedMalonylH.plastid.alphaWaterExchange=1;

T.retainedMalonylH.er.source='acetylStarterH';
T.retainedMalonylH.er.waterExchangeFraction=sourceVals(8);
T.retainedMalonylH.er.exchangeWaterSource='cytosolWater';
T.retainedMalonylH.er.alphaWaterExchange=1;

T.cer3AldehydeH.mode=P.cer3Source;

sourcePools=wax_hydrogen_source_pools_postmalonate_exactD(t,intracellular,T);
Rtraits=alphaVectorToTraitsLocal(reactionAlpha);
incorporatedPools=wax_apply_reaction_incorporation_exactD(sourcePools,Rtraits);
S=alkane_biosynthesis_h_stoichiometry_reaction_explicit( ...
    cfg.alkaneCarbonNumber,cfg.routeFractions);
W=wax_combine_hydrogen_classes_exactD(S,incorporatedPools);

rWax=1+W.deltaDMatrix/1000;
rPl=1+dPl/1000;
rCy=1+dCy/1000;

Bio=struct;
Bio.version='2026-08-29_parameterized_biochemical_bridge_v2';
Bio.timeDays=t;
Bio.parameterMode='explicitGSAParameters';
Bio.parameterInterpretation='explicit GSA/source/reaction parameter state';
Bio.sourceReferenceWater_permil=refWater;
Bio.sourceParameterNames=P.sourceParameterNames;
Bio.sourceParameterValues=sourceVals;
Bio.reactionParameterNames=P.reactionPriorDefinitions.parameterNames;
Bio.reactionAlphaValues=reactionAlpha;
Bio.sourceTraits=T;
Bio.sourcePools=sourcePools;
Bio.reactionTraits=Rtraits;
Bio.incorporatedPools=incorporatedPools;
Bio.stoichiometry=S;
Bio.alkaneCarbonNumber=cfg.alkaneCarbonNumber(:).';
Bio.wax=W;


if strcmpi(cfg.architectureVersion,'v2')
    Bio.waxStageExpanded=apply_wax_er_stage_states_exactD( ...
        W.deltaDMatrix,Bio.alkaneCarbonNumber,cfg.erStageStates);
else
    Tzero=wax_er_stage_state_traits();
    Tzero.stage6Delta_permil=0;
    Tzero.stage7Delta_permil=0;
    Bio.waxStageExpanded=apply_wax_er_stage_states_exactD( ...
        W.deltaDMatrix,Bio.alkaneCarbonNumber,Tzero);
end

Bio.architecture=struct;
Bio.architecture.version=cfg.architectureVersion;
Bio.architecture.baseWaxField='wax';
Bio.architecture.expandedWaxField='waxStageExpanded';
Bio.architecture.erStageStates=cfg.erStageStates;
Bio.architecture.zeroStateNestsV1=true;

Bio.apparentEpsilonReference.vsPlastidWater_permil=1000.*(rWax./rPl-1);
Bio.apparentEpsilonReference.vsCytosolWater_permil=1000.*(rWax./rCy-1);
Bio.apparentEpsilonReference.warning= ...
    ['Reference-only apparent epsilon; no unique scalar biochemical water ', ...
     'state exists when plastid and cytosol differ.'];
Bio.isCalibrated=false;
end

function E=prescribedEndmemberLocal(d)
E=struct('source','prescribed','deltaD_permil',d, ...
    'ratioNormalized',[],'alphaFormation',1);
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
