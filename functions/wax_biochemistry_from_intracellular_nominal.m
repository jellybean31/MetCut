function Bio = wax_biochemistry_from_intracellular_nominal( ...
    timeDays,intracellular,cfg)
%WAX_BIOCHEMISTRY_FROM_INTRACELLULAR_NOMINAL
% Deterministic bridge from the current exact-D intracellular-water state
% into the leaf-wax biosynthesis model.
%
%   plastid-derived biochemical source states
%       -> propagated from intracellular.deltaD.plastid
%
%   ER/cytosolic biochemical source states
%       -> propagated from intracellular.deltaD.cytosol
%
% The default parameter state is the numerical midpoint of the existing
% CORE biochemical DEVELOPMENT PRE-PRIORS. 
%
% Required production functions:
%   wax_integrated_biochemical_prior_definitions
%   wax_source_reference_offset_to_deltaD
%   wax_hydrogen_source_traits_postmalonate_null
%   wax_hydrogen_source_pools_postmalonate_exactD
%   wax_reaction_incorporation_traits_null
%   wax_apply_reaction_incorporation_exactD
%   alkane_biosynthesis_h_stoichiometry_reaction_explicit
%   wax_combine_hydrogen_classes_exactD

if nargin<3 || isempty(cfg), cfg=struct; end

if ~isfield(cfg,'parameterMode') || isempty(cfg.parameterMode)
    cfg.parameterMode='coreMidpointDiagnostic';
end
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

if ~strcmpi(cfg.parameterMode,'coreMidpointDiagnostic')
    error('Only coreMidpointDiagnostic is implemented in this first integrated bridge.');
end

requiredFunctions={ ...
    'wax_integrated_biochemical_prior_definitions', ...
    'wax_source_reference_offset_to_deltaD', ...
    'wax_hydrogen_source_traits_postmalonate_null', ...
    'wax_hydrogen_source_pools_postmalonate_exactD', ...
    'wax_reaction_incorporation_traits_null', ...
    'wax_apply_reaction_incorporation_exactD', ...
    'alkane_biosynthesis_h_stoichiometry_reaction_explicit', ...
    'wax_combine_hydrogen_classes_exactD'};

for k=1:numel(requiredFunctions)
    if exist(requiredFunctions{k},'file')~=2
        error('Required biochemical function is missing: %s.m',requiredFunctions{k});
    end
end

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

P=wax_integrated_biochemical_prior_definitions();

sourceMid=0.5.*(P.core.sourceLower+P.core.sourceUpper);
reactionMid=0.5.*( ...
    P.core.reaction.lowerAlpha(:).' + ...
    P.core.reaction.upperAlpha(:).');

% CER3 remains fixed at alpha=1 by the accepted prior definition.
reactionMid(7)=1;

refWater=P.sourceReferenceWater_permil;

% -------------------------------------------------------------------------
% Constant-alpha source propagation from compartment-appropriate waters.
% -------------------------------------------------------------------------

dAc= ...
    wax_source_reference_offset_to_deltaD( ...
        sourceMid(1),dPl,refWater);

dPNP= ...
    wax_source_reference_offset_to_deltaD( ...
        sourceMid(2),dPl,refWater);

dPNH= ...
    wax_source_reference_offset_to_deltaD( ...
        sourceMid(3),dPl,refWater);

dERNP= ...
    wax_source_reference_offset_to_deltaD( ...
        sourceMid(4),dCy,refWater);

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

T.plastidNADPH.waterExchangeFraction=sourceMid(5);
T.plastidNADPH.alphaWaterExchange=1;
T.plastidNADPH.exchangeWaterSource='plastidWater';

T.plastidNADH.waterExchangeFraction= ...
    P.plastidNADHWaterExchange.fixedFraction;

T.erNADPH.waterExchangeFraction=sourceMid(6);
T.erNADPH.alphaWaterExchange=1;
T.erNADPH.exchangeWaterSource='cytosolWater';

T.retainedMalonylH.plastid.source='acetylStarterH';
T.retainedMalonylH.plastid.waterExchangeFraction=sourceMid(7);
T.retainedMalonylH.plastid.exchangeWaterSource='plastidWater';
T.retainedMalonylH.plastid.alphaWaterExchange=1;

T.retainedMalonylH.er.source='acetylStarterH';
T.retainedMalonylH.er.waterExchangeFraction=sourceMid(8);
T.retainedMalonylH.er.exchangeWaterSource='cytosolWater';
T.retainedMalonylH.er.alphaWaterExchange=1;

T.cer3AldehydeH.mode=P.cer3Source;

sourcePools= ...
    wax_hydrogen_source_pools_postmalonate_exactD( ...
        t,intracellular,T);

Rtraits= ...
    alphaVectorToTraitsLocal(reactionMid);

incorporatedPools= ...
    wax_apply_reaction_incorporation_exactD( ...
        sourcePools,Rtraits);

S= ...
    alkane_biosynthesis_h_stoichiometry_reaction_explicit( ...
        cfg.alkaneCarbonNumber, ...
        cfg.routeFractions);

W= ...
    wax_combine_hydrogen_classes_exactD( ...
        S,incorporatedPools);

rWax=1+W.deltaDMatrix/1000;
rPl=1+dPl/1000;
rCy=1+dCy/1000;

epsVsPlastid=1000.*(rWax./rPl-1);
epsVsCytosol=1000.*(rWax./rCy-1);

Bio=struct;
Bio.timeDays=t;
Bio.parameterMode=cfg.parameterMode;
Bio.parameterInterpretation= ...
    'midpoint of core development pre-priors; diagnostic only; not calibrated';

Bio.sourceReferenceWater_permil=refWater;
Bio.sourcePropagationMode='constantAlphaFromReferenceOffset';
Bio.sourceParameterNames=P.sourceParameterNames;
Bio.sourceParameterValues=sourceMid;
Bio.reactionParameterNames=P.reactionPriorDefinitions.parameterNames;
Bio.reactionAlphaValues=reactionMid;

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

Bio.apparentEpsilonReference=struct;
Bio.apparentEpsilonReference.vsPlastidWater_permil=epsVsPlastid;
Bio.apparentEpsilonReference.vsCytosolWater_permil=epsVsCytosol;
Bio.apparentEpsilonReference.warning= ...
    ['Reference-only apparent epsilon. Wax receives H from multiple ', ...
     'compartment-specific biochemical source pools; no unique scalar ', ...
     'biosynthetic-water isotope state exists when plastid and cytosol differ.'];

Bio.isCalibrated=false;
Bio.modelRole='deterministic biochemical bridge for integrated budget development';
end

%% ========================================================================
function E=prescribedEndmemberLocal(d)
E=struct( ...
    'source','prescribed', ...
    'deltaD_permil',d, ...
    'ratioNormalized',[], ...
    'alphaFormation',1);
end

%% ========================================================================
function T=alphaVectorToTraitsLocal(X)

if numel(X)~=7
    error('Reaction alpha vector must contain seven values.');
end

T=wax_reaction_incorporation_traits_null();

T.plastid.KAR.alphaHydride=X(1);
T.plastid.ENR.alphaHydride=X(2);
T.plastid.ENR.alphaProton=X(3);

T.er.KCR.alphaHydride=X(4);
T.er.ECR.alphaHydride=X(5);
T.er.ECR.alphaProton=X(6);

T.cer3.alphaAldehydeHydride=X(7);
end
