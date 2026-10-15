function out = forced_surface_gb(uLeaf_m_s,characteristicDimension_m, ...
    momentumDiffusivity_m2_s,heatDiffusivity_m2_s,waterVaporDiffusivity_m2_s, ...
    Tair_K,pressure_Pa)
% FORCED_SURFACE_GB
%
% Generic production helper for forced-convection boundary-layer conductance.
%
% Surface order is strictly:
%   [ABAXIAL, ADAXIAL]
%
% REQUIRED INPUTS
% ---------------
% uLeaf_m_s                  1x2 or 2x1 leaf-level wind [m s^-1]
% characteristicDimension_m 1x2 or 2x1 characteristic dimension [m]
% momentumDiffusivity_m2_s  positive scalar [m^2 s^-1]
% heatDiffusivity_m2_s      positive scalar [m^2 s^-1]
% waterVaporDiffusivity_m2_s positive scalar [m^2 s^-1]
% Tair_K                     positive scalar [K]
% pressure_Pa                positive scalar [Pa]
%
% Re_s = u_s*d_s/Dm
%
% Nu_lam  = 0.6*Re^0.5
% Nu_turb = 0.032*Re^0.8
% Nu_forced = max(Nu_lam,Nu_turb)
%
% Sh_forced = Nu_forced*(Dh/Dw)^(1/3)
%
% gb_velocity = Dw*Sh_forced/d
%
% c_air = pressure/(R*Tair)
% gb_molar = gb_velocity*c_air

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
out.ShForced=[NaN NaN];
out.gbVelocity_m_s=[NaN NaN];
out.airMolarDensity_mol_m3=NaN;
out.gb_mol_m2_s=[NaN NaN];

% ---- type / shape guards ----
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

scalars=[momentumDiffusivity_m2_s,heatDiffusivity_m2_s, ...
    waterVaporDiffusivity_m2_s,Tair_K,pressure_Pa];

if ~all(isfinite([u,d,scalars]))
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
if ~(isscalar(waterVaporDiffusivity_m2_s) && waterVaporDiffusivity_m2_s>0)
    out.code="INVALID_VAPOR_DIFFUSIVITY";
    return
end
if ~(isscalar(Tair_K) && Tair_K>0)
    out.code="INVALID_AIR_TEMPERATURE";
    return
end
if ~(isscalar(pressure_Pa) && pressure_Pa>0)
    out.code="INVALID_PRESSURE";
    return
end

Re=u.*d./momentumDiffusivity_m2_s;
NuLam=0.6.*sqrt(Re);
NuTurb=0.032.*Re.^0.8;
NuForced=max(NuLam,NuTurb);

branch=repmat("LAMINAR_CORRELATION",1,2);
branch(NuTurb>NuLam)="TURBULENT_CORRELATION";
branch(Re==0)="ZERO_FORCED_CONVECTION";

Sh=NuForced.*(heatDiffusivity_m2_s/waterVaporDiffusivity_m2_s)^(1/3);
gbVel=waterVaporDiffusivity_m2_s.*Sh./d;
cAir=pressure_Pa/(Rgas*Tair_K);
gbMol=gbVel.*cAir;

if ~all(isfinite([Re,NuLam,NuTurb,NuForced,Sh,gbVel,cAir,gbMol]))
    out.code="NONFINITE_OUTPUT";
    return
end
if any(gbMol<0)
    out.code="NEGATIVE_OUTPUT";
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
out.ShForced=Sh;
out.gbVelocity_m_s=gbVel;
out.airMolarDensity_mol_m3=cAir;
out.gb_mol_m2_s=gbMol;
end
