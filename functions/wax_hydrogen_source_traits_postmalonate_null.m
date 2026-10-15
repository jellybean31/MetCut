function T = wax_hydrogen_source_traits_postmalonate_null()
%WAX_HYDROGEN_SOURCE_TRAITS_POSTMALONATE_NULL
% Null configuration for explicit acetyl / retained-malonyl / CER3 sources.

T = wax_hydrogen_source_traits_pathway_null();

legacyAcetyl = T.acetylMalonylH;

T.acetylStarterH = legacyAcetyl;

T.retainedMalonylH.plastid = legacyAcetyl;
T.retainedMalonylH.plastid.source = 'acetylStarterH';
T.retainedMalonylH.plastid.waterExchangeFraction = 0;
T.retainedMalonylH.plastid.exchangeWaterSource = 'plastidWater';
T.retainedMalonylH.plastid.alphaWaterExchange = 1;

T.retainedMalonylH.er = legacyAcetyl;
T.retainedMalonylH.er.source = 'acetylStarterH';
T.retainedMalonylH.er.waterExchangeFraction = 0;
T.retainedMalonylH.er.exchangeWaterSource = 'cytosolWater';
T.retainedMalonylH.er.alphaWaterExchange = 1;

T.cer3AldehydeH.mode = 'cytosolWater';
end
