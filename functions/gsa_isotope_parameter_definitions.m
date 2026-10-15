function D = gsa_isotope_parameter_definitions
%GSA_ISOTOPE_PARAMETER_DEFINITIONS
% GSA parameters supported by the current production model.
%
% SCREENING V4 / H4:
%   Samples 12 stable xylem/intracellular/biochemical dimensions.
%
% H4 change:
%   I08 is the grouped cytosol<->plastid exchange rate
%
%       k_CyP = Pm * (A_m/V_P)
%
% Temporarily FIXED at validated baseline because the spatial leaf ODE
% becomes prohibitively stiff across part of the empirical anatomy range:
%
%   I03 interveinalDistance_m        = 5e-4 m
%   I04 mesophyllThickness_m        = 1.5e-4 m
%   I05 effectiveFlowAreaFraction   = 0.02


D=struct( ...
    'id',{},'name',{},'low',{},'high',{},'distribution',{}, ...
    'stage',{},'field',{},'notes',{});

D(end+1)=p('I01','mobileTransitDays',0.5,5,'loguniform','xylem', ...
    'cfg.xylemIsotope.mobileTransitDays','days');

D(end+1)=p('I02','xylemStorageExchangeTauDays',2,20,'loguniform','xylem', ...
    'grouped lumen+living exchange tau','days');

D(end+1)=p('I07','cytosolWaterFraction',0.03,0.09,'uniform','intracellular', ...
    'cfg.intracellular.waterFractionOfSymplast.cytosol','fraction');

PX=intracellular_plastid_exchange_prior();

D(end+1)=p('I08','cytosolPlastidExchangeRate_s_inv', ...
    PX.groupedRate.low_s_inv,PX.groupedRate.high_s_inv, ...
    PX.groupedRate.distribution,'intracellular', ...
    'derived k = Pm*(A_m/V_P)', ...
    's^-1; grouped identifiable membrane-exchange rate');

D(end+1)=p('B01','acetylOffset_permil',-150,-50,'uniform','biochemistry', ...
    'sourceParameterValues(1)','per mil');

D(end+1)=p('B02','plastidNADPHOffset_permil',-300,-100,'uniform','biochemistry', ...
    'sourceParameterValues(2)','per mil');

D(end+1)=p('B03','plastidNADHOffset_permil',-200,-50,'uniform','biochemistry', ...
    'sourceParameterValues(3)','per mil');

D(end+1)=p('B04','erNADPHOffset_permil',-250,-50,'uniform','biochemistry', ...
    'sourceParameterValues(4)','per mil');

D(end+1)=p('B05','plastidNADPHExchangeFraction',0.4,1,'uniform','biochemistry', ...
    'sourceParameterValues(5)','fraction');

D(end+1)=p('B06','erNADPHExchangeFraction',0.4,1,'uniform','biochemistry', ...
    'sourceParameterValues(6)','fraction');

D(end+1)=p('B07','postMalonateExchangeFraction',0,1,'uniform','biochemistry', ...
    'sourceParameterValues(7:8), tied in first pass','fraction');

D(end+1)=p('B08','C16RouteFraction',0,0.5,'uniform','biochemistry', ...
    'cfg.biochemistry.routeFractions.C16; C18=1-C16','fraction');

end

function q=p(id,name,lo,hi,dist,stage,field,notes)
q=struct('id',id,'name',name,'low',lo,'high',hi, ...
    'distribution',dist,'stage',stage,'field',field,'notes',notes);
end
