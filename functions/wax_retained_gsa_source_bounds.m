
function G = wax_retained_gsa_source_bounds()
%WAX_RETAINED_GSA_SOURCE_BOUNDS
% Current retained primitive GSA biochemical parameter bounds.
%
% These are the retained GSA ranges used for the inverse/GSA design:
%   B01-B04 source offsets unchanged from core development envelope
%   B05 plastid NADPH-water exchange = 0.4-1
%   B06 ER/cytosolic NADPH-water exchange = 0.4-1
%   B07 grouped post-malonate exchange = 0-1

G=struct;

G.id={'B01','B02','B03','B04','B05','B06','B07'};
G.name={ ...
    'acetylStarterOffset_permil', ...
    'plastidNADPHOffset_permil', ...
    'plastidNADHOffset_permil', ...
    'erNADPHOffset_permil', ...
    'plastidNADPHWaterExchange', ...
    'erNADPHWaterExchange', ...
    'groupedPostMalonateExchange'};

G.lower=[-150 -200 -150 -175 0.4 0.4 0];
G.upper=[ -50  -50  -25  -25 1.0 1.0 1];

G.table=table( ...
    string(G.id(:)),string(G.name(:)),G.lower(:),G.upper(:), ...
    'VariableNames',{'ParameterID','Parameter','Lower','Upper'});

end
