function P = intracellular_plastid_exchange_prior
%INTRACELLULAR_PLASTID_EXCHANGE_PRIOR
% Physiological prior for cytosol<->plastid water exchange.
%
% The production model depends on membrane permeability and chloroplast
% geometry only through:
%
%   k_CyP = Pm * (A_m/V_P)
%
% Physiological screening range:
%   k_CyP = 10 .. 70 s^-1   (log-uniform)
%
% This range is relatively conservative:
%   - includes the current legacy baseline k=10 s^-1 (tau=0.1 s)
%   - is consistent with measured plant membrane permeabilities and
%     mesophyll chloroplast surface/volume ratios
%   - includes the order of magnitude implied by the classic in-vivo
%     chloroplast water residence-time measurement (~88 ms)
%
% For the production membranePhysical interface, k is represented using a
% fixed bookkeeping A/V reference and Pm=k/(A/V). This fixed A/V is NOT an
% independently inferred anatomical parameter.

P=struct;
P.version='2026-08-29_intracellular_plastid_exchange_prior_H4_v1';

P.groupedRate=struct;
P.groupedRate.name='cytosolPlastidExchangeRate';
P.groupedRate.low_s_inv=10;
P.groupedRate.high_s_inv=70;
P.groupedRate.distribution='loguniform';
P.groupedRate.baseline_s_inv=10;

P.referenceGeometry=struct;
P.referenceGeometry.areaToVolume_m_inv=3e6;
P.referenceGeometry.status= ...
    'bookkeeping reference for mapping grouped k into membranePhysical mode';
P.referenceGeometry.notSeparatelyInferred=true;

P.derived=struct;
P.derived.tauLow_s=1/P.groupedRate.high_s_inv;
P.derived.tauHigh_s=1/P.groupedRate.low_s_inv;
P.derived.PmLow_m_s=P.groupedRate.low_s_inv / ...
    P.referenceGeometry.areaToVolume_m_inv;
P.derived.PmHigh_m_s=P.groupedRate.high_s_inv / ...
    P.referenceGeometry.areaToVolume_m_inv;

P.registryContext=struct;
P.registryContext.PmScreenLow_m_s=5e-6;
P.registryContext.PmScreenHigh_m_s=2e-5;
P.registryContext.note= ...
    ['Registry Pm range retained as physiological context; grouped k is ', ...
     'sampled because Pm and A/V are structurally non-identifiable separately.'];

P.provenance=struct;
P.provenance.chloroplastWaterPermeability= ...
    'McCain & Markley 1985, FEBS Lett, DOI 10.1016/0014-5793(85)80809-7';
P.provenance.plantProtoplastPermeability= ...
    'Sommer et al. 2007, J Membr Biol, DOI 10.1007/s00232-007-9011-6';
P.provenance.riceChloroplastGeometry= ...
    'Oikawa et al. 2017, Ann Bot, chloroplast S/V 3.43 um^-1';
P.provenance.C4ChloroplastGeometry= ...
    'Lee et al. 2023, New Phytol, mesophyll chloroplast S/V ~2.4-2.6 um^-1';

P.interpretation= ...
    ['Use k_CyP as the first-pass inferable/sensitivity parameter. ', ...
     'Map back to Pm only when species-specific chloroplast A/V is known.'];

end
