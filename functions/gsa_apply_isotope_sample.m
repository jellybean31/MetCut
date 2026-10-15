function S = gsa_apply_isotope_sample(Mreference,sample)
%GSA_APPLY_ISOTOPE_SAMPLE Build replay overrides from one GSA sample row.
%
% H4:
%   I08 samples grouped cytosol<->plastid rate
%
%       k_CyP = Pm*(A_m/V_P)
%
%   and maps it into the production membranePhysical interface using the
%   fixed bookkeeping A/V reference returned by
%   intracellular_plastid_exchange_prior.
%
% sample may be a one-row table or scalar struct created from
% gsa_sample_lhs. Returns a normalized override structure.

if istable(sample)
    if height(sample)~=1, error('sample table must contain one row.'); end
    sample=table2struct(sample);
end
if ~isstruct(sample) || ~isscalar(sample)
    error('sample must be a scalar struct or one-row table.');
end

cfg=Mreference.config;
P=wax_integrated_biochemical_prior_definitions();
sourceVals=0.5.*(P.core.sourceLower+P.core.sourceUpper);
reactionAlpha=0.5.*(P.core.reaction.lowerAlpha(:).' + P.core.reaction.upperAlpha(:).');
reactionAlpha(7)=1;

S=struct;
S.cfg=cfg;
S.sourceVals=sourceVals;
S.reactionAlpha=reactionAlpha;
S.changedStages={};
S.derived=struct;

getv=@(id) getByPrefixLocal(sample,id);

v=getv('I01');
if ~isempty(v)
    S.cfg.xylemIsotope.mobileTransitDays=v;
    S.changedStages{end+1}='xylem';
end

v=getv('I02');
if ~isempty(v)
    S.cfg.xylemIsotope.storageLumenExchangeTauDays=v;
    S.cfg.xylemIsotope.storageLivingExchangeTauDays=v;
    S.changedStages{end+1}='xylem';
end

v=getv('I03');
if ~isempty(v)
    S.cfg.leafIsotope.anatomy.interveinalDistance_m=v;
    S.changedStages{end+1}='leaf';
end

v=getv('I04');
if ~isempty(v)
    S.cfg.leafIsotope.anatomy.mesophyllThickness_m=v;
    S.changedStages{end+1}='leaf';
end

v=getv('I05');
if ~isempty(v)
    S.cfg.leafIsotope.anatomy.effectiveFlowAreaFraction=v;
    S.changedStages{end+1}='leaf';
end

v=getv('I07');
if ~isempty(v)
    S.cfg.intracellular.waterFractionOfSymplast.cytosol=v;
    S.changedStages{end+1}='intracellular';
end

v=getv('I08');
if ~isempty(v)

    PX=intracellular_plastid_exchange_prior();
    AV=PX.referenceGeometry.areaToVolume_m_inv;

    if ~isscalar(v) || ~isfinite(v) || v<=0
        error('I08 grouped cytosol-plastid exchange rate must be positive.');
    end

    S.cfg.intracellular.cytosolPlastidExchange.mode='membranePhysical';
    S.cfg.intracellular.cytosolPlastidExchange. ...
        membraneAreaToPlastidWaterVolume_m_inv=AV;
    S.cfg.intracellular.cytosolPlastidExchange. ...
        membraneWaterPermeability_m_s=v/AV;

    S.derived.I08_kCyP_s_inv=v;
    S.derived.I08_tauCyP_s=1/v;
    S.derived.I08_referenceAreaToVolume_m_inv=AV;
    S.derived.I08_impliedPm_m_s=v/AV;

    S.changedStages{end+1}='intracellular';
end

map={'B01',1;'B02',2;'B03',3;'B04',4;'B05',5;'B06',6};
for k=1:size(map,1)
    v=getv(map{k,1});
    if ~isempty(v)
        S.sourceVals(map{k,2})=v;
        S.changedStages{end+1}='biochemistry';
    end
end

v=getv('B07');
if ~isempty(v)
    S.sourceVals(7)=v;
    S.sourceVals(8)=v;
    S.changedStages{end+1}='biochemistry';
end

v=getv('B08');
if ~isempty(v)
    if v<0 || v>1, error('C16 route fraction must lie in [0,1].'); end
    S.cfg.biochemistry.routeFractions.C16=v;
    S.cfg.biochemistry.routeFractions.C18=1-v;
    S.changedStages{end+1}='biochemistry';
end

S.changedStages=unique(S.changedStages,'stable');
end

function v=getByPrefixLocal(s,id)
v=[];
f=fieldnames(s);
for k=1:numel(f)
    if startsWith(f{k},[id '_']) || strcmp(f{k},id)
        v=s.(f{k});
        return
    end
end
end
