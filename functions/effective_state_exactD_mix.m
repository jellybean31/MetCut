
function d = effective_state_exactD_mix(dSource,dWater,fExchange)
%EFFECTIVE_STATE_EXACTD_MIX Exact D-atom-fraction mixing for H exchange.
%
% dSource and dWater are deltaD in per mil.
% fExchange is in [0,1].

if any(fExchange(:)<0 | fExchange(:)>1)
    error('fExchange must lie in [0,1].');
end

Rstd=155.76e-6; % VSMOW D/H used by the wax source-pool functions

rS=1+dSource./1000;
rW=1+dWater./1000;

RS=Rstd.*rS;
RW=Rstd.*rW;

xS=RS./(1+RS);
xW=RW./(1+RW);

x=(1-fExchange).*xS+fExchange.*xW;

R=x./(1-x);
r=R./Rstd;

d=1000.*(r-1);

end
