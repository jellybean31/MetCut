function out = forced_surface_heat( ...
    uLeaf_m_s,characteristicDimension_m,momentumDiffusivity_m2_s, ...
    heatDiffusivity_m2_s,Tair_K,Tleaf_K,pressure_Pa)
% FORCED_SURFACE_HEAT
%
% Generic production helper for forced-convection sensible heat.
%
% Surface order:
%   [ABAXIAL, ADAXIAL]
%
% REQUIRED INPUTS
% ---------------
% uLeaf_m_s                  1x2 or 2x1 leaf-level wind [m s^-1]
% characteristicDimension_m 1x2 or 2x1 characteristic dimension [m]
% momentumDiffusivity_m2_s  positive scalar [m^2 s^-1]
% heatDiffusivity_m2_s      positive scalar [m^2 s^-1]
% Tair_K                     positive scalar [K]
% Tleaf_K                    positive scalar [K]
% pressure_Pa                positive scalar [Pa]
%
% -------------
% Re = u*d/Dm
% Nu_lam = 0.6*Re^0.5
% Nu_turb = 0.032*Re^0.8
% Nu_forced = max(Nu_lam,Nu_turb)
%
% Tfilm = 0.5*(Tleaf+Tair)
% gH_velocity = Dh*Nu_forced/d
% c_air = P/(R*Tfilm)
% gH_molar = gH_velocity*c_air
% cp_air_molar = 29.07 + 0.003*(Tfilm-273.15)
% H_surface = cp_air_molar*gH_molar*(Tleaf-Tair)
%
% SIGN
% ----
% H_surface > 0 means heat leaves the leaf.


Rgas=8.3144598;

out=struct();
out.success=false;
out.code="UNSET";
out.surfaceOrder=["ABAXIAL","ADAXIAL"];
out.uLeaf_m_s=[NaN NaN];
out.characteristicDimension_m=[NaN NaN];
out.Re=[NaN NaN];
out.NuLaminar=[NaN NaN];
out.NuTurbulent=[NaN NaN];
out.NuForced=[NaN NaN];
out.selectedBranch=["",""];
out.Tair_K=NaN;
out.Tleaf_K=NaN;
out.Tfilm_K=NaN;
out.cpAirMolar_J_mol_K=NaN;
out.gHVelocity_m_s=[NaN NaN];
out.airMolarDensity_mol_m3=NaN;
out.gHMolar_mol_m2_s=[NaN NaN];
out.Hsurface_W_m2=[NaN NaN];

if ~(isnumeric(uLeaf_m_s) && isreal(uLeaf_m_s) && numel(uLeaf_m_s)==2)
    out.code="INVALID_WIND_SHAPE";
    return
end
if ~(isnumeric(characteristicDimension_m) && isreal(characteristicDimension_m) && ...
        numel(characteristicDimension_m)==2)
    out.code="INVALID_DIMENSION_SHAPE";
    return
end

u=reshape(double(uLeaf_m_s),1,2);
d=reshape(double(characteristicDimension_m),1,2);

vals=[u,d,momentumDiffusivity_m2_s,heatDiffusivity_m2_s, ...
    Tair_K,Tleaf_K,pressure_Pa];

if ~all(isfinite(vals))
    out.code="NONFINITE_INPUT";
    return
end
if any(u<0)
    out.code="NEGATIVE_LEAF_WIND";
    return
end
if any(d<=0)
    out.code="NONPOSITIVE_CHARACTERISTIC_DIMENSION";
    return
end
if ~(isscalar(momentumDiffusivity_m2_s) && momentumDiffusivity_m2_s>0)
    out.code="INVALID_MOMENTUM_DIFFUSIVITY";
    return
end
if ~(isscalar(heatDiffusivity_m2_s) && heatDiffusivity_m2_s>0)
    out.code="INVALID_HEAT_DIFFUSIVITY";
    return
end
if ~(isscalar(Tair_K) && Tair_K>0)
    out.code="INVALID_AIR_TEMPERATURE";
    return
end
if ~(isscalar(Tleaf_K) && Tleaf_K>0)
    out.code="INVALID_LEAF_TEMPERATURE";
    return
end
if ~(isscalar(pressure_Pa) && pressure_Pa>0)
    out.code="INVALID_PRESSURE";
    return
end

Tfilm=0.5*(Tair_K+Tleaf_K);
TfilmC=Tfilm-273.15;
cpAir=29.07+0.003*TfilmC;

if ~(isfinite(cpAir) && cpAir>0)
    out.code="INVALID_AIR_HEAT_CAPACITY";
    return
end

Re=u.*d./momentumDiffusivity_m2_s;
NuLam=0.6.*sqrt(Re);
NuTurb=0.032.*Re.^0.8;
NuForced=max(NuLam,NuTurb);

branch=repmat("LAMINAR_CORRELATION",1,2);
branch(NuTurb>NuLam)="TURBULENT_CORRELATION";
branch(Re==0)="ZERO_FORCED_CONVECTION";

gHVel=heatDiffusivity_m2_s.*NuForced./d;
cAir=pressure_Pa/(Rgas*Tfilm);
gHMol=gHVel.*cAir;
H=cpAir.*gHMol.*(Tleaf_K-Tair_K);

if ~all(isfinite([Re,NuLam,NuTurb,NuForced,gHVel,cAir,gHMol,H]))
    out.code="NONFINITE_OUTPUT";
    return
end
if any(gHVel<0) || any(gHMol<0)
    out.code="NEGATIVE_CONDUCTANCE_OUTPUT";
    return
end

out.success=true;
out.code="OK";
out.uLeaf_m_s=u;
out.characteristicDimension_m=d;
out.Re=Re;
out.NuLaminar=NuLam;
out.NuTurbulent=NuTurb;
out.NuForced=NuForced;
out.selectedBranch=branch;
out.Tair_K=Tair_K;
out.Tleaf_K=Tleaf_K;
out.Tfilm_K=Tfilm;
out.cpAirMolar_J_mol_K=cpAir;
out.gHVelocity_m_s=gHVel;
out.airMolarDensity_mol_m3=cAir;
out.gHMolar_mol_m2_s=gHMol;
out.Hsurface_W_m2=H;
end
