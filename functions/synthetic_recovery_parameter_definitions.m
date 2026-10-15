
function D = synthetic_recovery_parameter_definitions
%SYNTHETIC_RECOVERY_PARAMETER_DEFINITIONS
% Parameter definitions for staged synthetic recovery experiments.
%
% Transport:
%   I01 mobileTransitDays
%   I02 xylemStorageExchangeTauDays
%
% Biochemistry:
%   B01 acetylOffset_permil
%   B02 plastidNADPHOffset_permil
%   B03 plastidNADHOffset_permil
%   B04 erNADPHOffset_permil
%   B05 plastidNADPHExchangeFraction
%   B06 erNADPHExchangeFraction
%   B07 postMalonateExchangeFraction
%

Dall=gsa_isotope_parameter_definitions();
ids=string({Dall.id});

keepTransport=["I01","I02"];
keepBio=["B01","B02","B03","B04","B05","B06","B07"];

D=struct;
D.transport=orderedSubsetLocal(Dall,ids,keepTransport);
D.biochemistry=orderedSubsetLocal(Dall,ids,keepBio);
D.integrated=[D.transport D.biochemistry];

end

function Dout=orderedSubsetLocal(Dall,ids,keep)
Dout=repmat(Dall(1),1,numel(keep));
for k=1:numel(keep)
    j=find(ids==keep(k),1);
    if isempty(j)
        error('Parameter %s not found in GSA definitions.',keep(k));
    end
    Dout(k)=Dall(j);
end
end
