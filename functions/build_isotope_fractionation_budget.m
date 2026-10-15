function B = build_isotope_fractionation_budget(M,deltaDPrecipInterval_permil)
%BUILD_ISOTOPE_FRACTIONATION_BUDGET
% Assemble an end-to-end isotope-state and apparent-transfer budget from
% the integrated model object.
%
% CRITICAL INTERPRETATION
% -----------------------
% "apparentTransfer_permil" is:
%
%   1000 * (R_to/R_from - 1)
%
% Source weighting, mixing, storage memory, transport lag and evaporative
% enrichment can all change isotope states.

t=M.timeDays(:);
nT=numel(t);

S=M.soilIsotopes;
R=M.rootBoundary;
X=M.xylemIsotopes;
L=M.leafIsotopes;
I=M.intracellular;
Bio=M.biochemistry;

%% ========================================================================
% State definitions
% ========================================================================

states=struct;

% Precipitation is an interval forcing. Convert to ending-node form only
% for plotting/state bookkeeping; node 1 uses first interval value.
dP=deltaDPrecipInterval_permil(:);
if numel(dP)~=nT-1
    error('deltaDPrecipInterval_permil must have nTime-1 values.');
end
states.precipitation.deltaD_permil=[dP(1);dP];
states.precipitation.activeMask=[false; ...
    S.transfer.infiltration_mol_m2(:)>0];
states.precipitation.definition= ...
    'Prescribed precipitation isotope forcing; only active where actual infiltration > 0.';

% Exact inventory-weighted bulk resolved-soil isotope state.
Wsoil=S.waterAmount_mol_m2;
xSoil=S.deuteriumAtomFraction;
xBulk=sum(Wsoil.*xSoil,1)./sum(Wsoil,1);
states.soilBulk.deltaD_permil=atomFractionToDeltaLocal(xBulk,M.config.Rstandard_DH).';
states.soilBulk.definition= ...
    'Exact water-inventory-weighted isotope composition across all resolved soil cells.';

states.rootSource.deltaD_permil=R.deltaDNode_permil(:);
states.rootSource.definition=R.definition;

states.xylemDelivered.deltaD_permil=X.leafBoundary.deltaDWater(:);
states.xylemDelivered.definition= ...
    'Exact-D xylem water delivered to the leaf boundary after mobile transit and storage exchange.';

states.leafApoplast.deltaD_permil=L.deltaDApoplast(:);
states.leafApoplast.definition= ...
    'Finite leaf apoplastic reservoir isotope composition.';

states.leafSymplastBulk.deltaD_permil=L.deltaDSymplastBulk(:);
states.leafSymplastBulk.definition= ...
    'Water-inventory-weighted bulk spatial symplast/mesophyll isotope composition.';

states.evaporationSite.deltaD_permil=L.deltaDEvaporationSite(:);
states.evaporationSite.definition= ...
    'Massless spatial evaporation-site isotope state solved with atmospheric fractionating exchange.';

states.cytosol.deltaD_permil=I.deltaD.cytosol(:);
states.cytosol.definition= ...
    'Exact-D intracellular cytosol tracer-water state; no intrinsic fractionation in intracellular exchange.';

states.plastid.deltaD_permil=I.deltaD.plastid(:);
states.plastid.definition= ...
    'Exact-D intracellular plastid tracer-water state; no intrinsic fractionation in intracellular exchange.';

chains=Bio.alkaneCarbonNumber(:);
waxDelta=Bio.wax.deltaDMatrix;

for h=1:numel(chains)
    nm=sprintf('waxC%d',chains(h));
    states.(nm).deltaD_permil=waxDelta(h,:).';
    states.(nm).definition=sprintf( ...
        'Reaction D inventory C%d n-alkane product.',chains(h));
end


% Exact stoichiometric C29 mixture BEFORE reaction-specific incorporation.
% This is the correct intermediate for separating:
%
%   intracellular water
%       -> biochemical source state / exchange
%       -> reaction-specific incorporation
%       -> final wax
%
% It remains a mixed plastid + cytosol/ER source state, so both intracellular
% reference waters are retained separately in the budget.
idx29=find(chains==29,1);
if isempty(idx29)
    error('C29 is required for the biochemical budget decomposition.');
end

sourceNames29=Bio.wax.sourceClassNames;
frac29=Bio.wax.compounds(idx29).sourceFractionVector(:);
preClasses29=Bio.incorporatedPools.preIncorporationClasses;

xPre29=zeros(numel(sourceNames29),nT);
for j=1:numel(sourceNames29)
    nm=sourceNames29{j};
    xPre29(j,:)=preClasses29.(nm).DAtomFraction(:).';
end

xPreMix29=frac29.'*xPre29;
dPreMix29=atomFractionToDeltaLocal(xPreMix29,Bio.wax.R_VSMOW).';

states.preReactionC29Mix.deltaD_permil=dPreMix29;
states.preReactionC29Mix.definition= ...
    ['Exact stoichiometric C29 H-source mixture after source-state formation ', ...
     'and exchange but before reaction-specific incorporation alphas.'];

%% ========================================================================
% Transition metadata
% ========================================================================

from={ ...
    'precipitation'; ...
    'soilBulk'; ...
    'rootSource'; ...
    'xylemDelivered'; ...
    'leafApoplast'; ...
    'leafSymplastBulk'; ...
    'leafSymplastBulk'; ...
    'cytosol'; ...
    'plastid'; ...
    'cytosol'; ...
    'preReactionC29Mix'};

to={ ...
    'soilBulk'; ...
    'rootSource'; ...
    'xylemDelivered'; ...
    'leafApoplast'; ...
    'leafSymplastBulk'; ...
    'evaporationSite'; ...
    'cytosol'; ...
    'plastid'; ...
    'preReactionC29Mix'; ...
    'preReactionC29Mix'; ...
    'waxC29'};

mechanism={ ...
    'infiltration + soil mixing/storage + vertical routing'; ...
    'root-access/source weighting + common-root mixing'; ...
    'mobile xylem advection + storage/recharge/release + residence-time memory'; ...
    'hydraulic delivery into finite apoplastic storage'; ...
    'coupled leaf evaporative-enrichment field: spatial transport/storage responding to the fractionating evaporation boundary'; ...
    'residual spatial gradient within the same coupled leaf evaporative-enrichment field'; ...
    'gross nonfractionating tracer exchange + storage adjustment'; ...
    'gross nonfractionating tracer exchange + storage adjustment'; ...
    'biochemical source-state formation/exchange; plastid water is a reference'; ...
    'biochemical source-state formation/exchange; cytosol water is a reference'; ...
    'reaction-specific biochemical incorporation from the pre-reaction C29 source mixture'};

intrinsicAlpha={ ...
    'not represented by one alpha'; ...
    'root uptake intrinsic alpha = 1 in current model'; ...
    'no intrinsic xylem fractionation unless isotopeSelectivity is explicitly enabled'; ...
    'no intrinsic fractionation assigned'; ...
    'not independently separable from the coupled fractionating evaporation boundary'; ...
    'same coupled fractionating evaporation physics; not the total leaf fractionation'; ...
    'alpha = 1'; ...
    'alpha = 1'; ...
    'multiple source-state and exchange alphas'; ...
    'multiple source-state and exchange alphas'; ...
    'multiple reaction-specific incorporation alphas'};

transitionMetadata=table( ...
    from,to,mechanism,intrinsicAlpha, ...
    'VariableNames',{ ...
    'FromState','ToState','Mechanism','IntrinsicFractionationStatus'});

%% ========================================================================
% Apparent transfer arrays
% ========================================================================

transfer=struct;

for k=1:height(transitionMetadata)

    a=transitionMetadata.FromState{k};
    b=transitionMetadata.ToState{k};

    dA=states.(a).deltaD_permil(:);
    dB=states.(b).deltaD_permil(:);

    if numel(dA)~=nT || numel(dB)~=nT
        error('Budget state length mismatch for %s -> %s.',a,b);
    end

    rA=1+dA/1000;
    rB=1+dB/1000;

    x=1000.*(rB./rA-1);

    fieldName=matlab.lang.makeValidName([a '_to_' b]);

    transfer.(fieldName)=struct;
    transfer.(fieldName).apparentTransfer_permil=x;
    transfer.(fieldName).deltaDifference_permil=dB-dA;
    transfer.(fieldName).fromDeltaD_permil=dA;
    transfer.(fieldName).toDeltaD_permil=dB;
    transfer.(fieldName).mechanism=mechanism{k};
    transfer.(fieldName).intrinsicFractionationStatus=intrinsicAlpha{k};
end

% Precipitation same-time transfer is only meaningful while infiltration is
% active. Explicitly mask it outside those nodes.
f=matlab.lang.makeValidName('precipitation_to_soilBulk');
transfer.(f).apparentTransfer_permil( ...
    ~states.precipitation.activeMask)=NaN;
transfer.(f).deltaDifference_permil( ...
    ~states.precipitation.activeMask)=NaN;

%% ========================================================================
% Closure / diagnostic summary
% ========================================================================

closure=struct;

closure.maxSoilDResidual_molD_m2= ...
    max(abs(S.stepIsotopeBalanceResidual_mol_D_m2),[],'omitnan');

closure.maxXylemDResidual_molD_m2= ...
    max(abs(X.deuteriumInventoryClosureResidual_molD_m2),[],'omitnan');

closure.maxLeafDResidual_molD_m2= ...
    max(abs(L.DInventoryClosureResidual_molD_m2),[],'omitnan');

closure.maxIntracellularDResidual_molD_m2= ...
    max(abs(I.wholeSystemDClosureResidual_molD_m2),[],'omitnan');

%% ========================================================================
% Biochemical source-class snapshot
% ========================================================================

biochemical=struct;

biochemical.sourceClassNames=Bio.sourcePools.sourceClassNames;
biochemical.referenceWarning=Bio.apparentEpsilonReference.warning;
biochemical.C29_vs_plastidWater_apparentEpsilon_permil= ...
    Bio.apparentEpsilonReference.vsPlastidWater_permil(2,:).';
biochemical.C29_vs_cytosolWater_apparentEpsilon_permil= ...
    Bio.apparentEpsilonReference.vsCytosolWater_permil(2,:).';

for k=1:numel(Bio.sourcePools.sourceClassNames)
    nm=Bio.sourcePools.sourceClassNames{k};
    biochemical.sourceDeltaD_permil.(nm)= ...
        Bio.sourcePools.classes.(nm).deltaD(:);
end

% Exact biochemical decomposition.
dPl=states.plastid.deltaD_permil(:);
dCy=states.cytosol.deltaD_permil(:);
dPre=states.preReactionC29Mix.deltaD_permil(:);
dWax=states.waxC29.deltaD_permil(:);

biochemical.preReactionC29Mix_deltaD_permil=dPre;

biochemical.C29_sourceStateEffect_vsPlastid_permil= ...
    apparentTransferLocal(dPl,dPre);
biochemical.C29_sourceStateEffect_vsCytosol_permil= ...
    apparentTransferLocal(dCy,dPre);

biochemical.C29_reactionIncorporationEffect_permil= ...
    apparentTransferLocal(dPre,dWax);

biochemical.C29_totalEffect_vsPlastid_permil= ...
    apparentTransferLocal(dPl,dWax);
biochemical.C29_totalEffect_vsCytosol_permil= ...
    apparentTransferLocal(dCy,dWax);

reconstructedPlastid=1000.*( ...
    (1+biochemical.C29_sourceStateEffect_vsPlastid_permil/1000).* ...
    (1+biochemical.C29_reactionIncorporationEffect_permil/1000)-1);

reconstructedCytosol=1000.*( ...
    (1+biochemical.C29_sourceStateEffect_vsCytosol_permil/1000).* ...
    (1+biochemical.C29_reactionIncorporationEffect_permil/1000)-1);

biochemical.C29_multiplicativeClosureResidual_vsPlastid_permil= ...
    reconstructedPlastid-biochemical.C29_totalEffect_vsPlastid_permil;

biochemical.C29_multiplicativeClosureResidual_vsCytosol_permil= ...
    reconstructedCytosol-biochemical.C29_totalEffect_vsCytosol_permil;

%% ========================================================================
% Output
% ========================================================================

B=struct;
B.version='2026-08-19_isotope_fractionation_budget_v2';
B.timeDays=t;
B.states=states;
B.transitionMetadata=transitionMetadata;
B.transfer=transfer;
B.biochemical=biochemical;
B.closure=closure;

B.interpretiveRules={ ...
    'Use apparentTransfer_permil as a state-to-state diagnostic.', ...
    'Do not add per-mil transitions to obtain rain-to-wax fractionation.', ...
    'Storage/memory processes require causal/lag attribution.', ...
    'Root uptake is nonfractionating in the current model; isotope change is source weighting.', ...
    'Intracellular exchange is nonfractionating; isotope change is residence-time/storage memory.', ...
    'Wax has multiple plastid/cytosol biochemical H sources, so no unique scalar biosynthetic-water epsilon exists once those pools diverge.', ...
    'Separate biochemical source-state/exchange effects from reaction incorporation using the exact preReactionC29Mix intermediate.', ...
    'Soil evaporative isotope fractionation is not yet connected and must be added before final end-to-end validation.'};

end


%% ========================================================================
function x=apparentTransferLocal(dFrom,dTo)
x=1000.*((1+dTo(:)/1000)./(1+dFrom(:)/1000)-1);
end

%% ========================================================================
function delta=atomFractionToDeltaLocal(xD,Rstd)
R=xD./(1-xD);
delta=1000.*(R./Rstd-1);
end
