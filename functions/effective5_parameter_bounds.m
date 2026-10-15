
function D = effective5_parameter_bounds(intracellular)
%EFFECTIVE5_PARAMETER_BOUNDS
% Effective-state bounds derived from GSA primitive ranges.


P=wax_integrated_biochemical_prior_definitions();
G=wax_retained_gsa_source_bounds();

lo=G.lower;
hi=G.upper;

dPl=intracellular.deltaD.plastid(end);
dCy=intracellular.deltaD.cytosol(end);

starterVals=[];
retVals=[];
pnpVals=[];
pnhVals=[];
erVals=[];

% E01: acetyl starter, B01.
for b=[lo(1) hi(1)]
    dPre=wax_source_reference_offset_to_deltaD( ...
        b,dPl,P.sourceReferenceWater_permil);
    starterVals(end+1)=dPre; 

    % E02: retained post-malonate precursor, grouped B07.
    for f=[lo(7) hi(7)]
        retVals(end+1)=effective_state_exactD_mix(dPre,dPl,f); 
    end
end

% E03: plastid NADPH, B02 + B05.
for b=[lo(2) hi(2)]
    dPre=wax_source_reference_offset_to_deltaD( ...
        b,dPl,P.sourceReferenceWater_permil);
    for f=[lo(5) hi(5)]
        pnpVals(end+1)=effective_state_exactD_mix(dPre,dPl,f);
    end
end

% E04: plastid NADH, B03; exchange remains fixed at production setting.
for b=[lo(3) hi(3)]
    dPre=wax_source_reference_offset_to_deltaD( ...
        b,dPl,P.sourceReferenceWater_permil);
    f=P.plastidNADHWaterExchange.fixedFraction;
    pnhVals(end+1)=effective_state_exactD_mix(dPre,dPl,f); 
end

% E05: ER/cytosolic NADPH, B04 + B06.
for b=[lo(4) hi(4)]
    dPre=wax_source_reference_offset_to_deltaD( ...
        b,dCy,P.sourceReferenceWater_permil);
    for f=[lo(6) hi(6)]
        erVals(end+1)=effective_state_exactD_mix(dPre,dCy,f); 
    end
end

D=struct( ...
    'id',{'E01','E02','E03','E04','E05'}, ...
    'name',{'acetylStarter_permil', ...
            'retainedPrecursor_permil', ...
            'plastidNADPH_permil', ...
            'plastidNADH_permil', ...
            'erNADPH_permil'}, ...
    'low',{min(starterVals),min(retVals),min(pnpVals),min(pnhVals),min(erVals)}, ...
    'high',{max(starterVals),max(retVals),max(pnpVals),max(pnhVals),max(erVals)});

end
