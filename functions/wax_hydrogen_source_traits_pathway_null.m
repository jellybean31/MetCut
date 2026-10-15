function T = wax_hydrogen_source_traits_pathway_null()
%WAX_HYDROGEN_SOURCE_TRAITS_PATHWAY_NULL
%
% Complete pathway-explicit NULL configuration for regression testing.
%
% Every biochemical H end-member is made isotopically identical to its
% local intracellular water, all formation alpha values are 1, and active-H
% waterExchangeFraction is 0.

T = struct;

%% Acetyl/malonyl inherited H
T.acetylMalonylH = endmember( ...
    'plastidWater',[],[],1);

%% Plastid NADPH
T.plastidNADPH.pathwayFractions = struct( ...
    'FNR',1, ...
    'OPPP',0, ...
    'NADPME',0);

T.plastidNADPH.pathways.FNR = endmember( ...
    'plastidWater',[],[],1);

T.plastidNADPH.pathways.OPPP.G6PDH = endmember( ...
    'plastidWater',[],[],1);

T.plastidNADPH.pathways.OPPP.SixPGDH = endmember( ...
    'plastidWater',[],[],1);

T.plastidNADPH.pathways.NADPME = endmember( ...
    'plastidWater',[],[],1);

T.plastidNADPH.waterExchangeFraction = 0;
T.plastidNADPH.exchangeWaterSource = 'plastidWater';
T.plastidNADPH.alphaWaterExchange = 1;

%% Plastid NADH
T.plastidNADH.pathwayFractions = struct( ...
    'pPDH',1, ...
    'GAPDH',0, ...
    'MDH',0, ...
    'PGDH',0);

T.plastidNADH.pathways.pPDH = endmember( ...
    'plastidWater',[],[],1);

T.plastidNADH.pathways.GAPDH = endmember( ...
    'plastidWater',[],[],1);

T.plastidNADH.pathways.MDH = endmember( ...
    'plastidWater',[],[],1);

T.plastidNADH.pathways.PGDH = endmember( ...
    'plastidWater',[],[],1);

T.plastidNADH.waterExchangeFraction = 0;
T.plastidNADH.exchangeWaterSource = 'plastidWater';
T.plastidNADH.alphaWaterExchange = 1;

%% ER/cytosolic NADPH
T.erNADPH.pathwayFractions = struct( ...
    'OPPP',1, ...
    'ICDH',0, ...
    'npGAPDH',0, ...
    'NADPME',0);

T.erNADPH.pathways.OPPP.G6PDH = endmember( ...
    'cytosolWater',[],[],1);

T.erNADPH.pathways.OPPP.SixPGDH = endmember( ...
    'cytosolWater',[],[],1);

T.erNADPH.pathways.ICDH = endmember( ...
    'cytosolWater',[],[],1);

T.erNADPH.pathways.npGAPDH = endmember( ...
    'cytosolWater',[],[],1);

T.erNADPH.pathways.NADPME = endmember( ...
    'cytosolWater',[],[],1);

T.erNADPH.waterExchangeFraction = 0;
T.erNADPH.exchangeWaterSource = 'cytosolWater';
T.erNADPH.alphaWaterExchange = 1;

%% Terminal H
T.terminalH.mode = 'cytosolWater';

end

%% ------------------------------------------------------------------------
function E = endmember(source,deltaD,ratio,alpha)

E = struct;
E.source = source;
E.deltaD_permil = deltaD;
E.ratioNormalized = ratio;
E.alphaFormation = alpha;

end
