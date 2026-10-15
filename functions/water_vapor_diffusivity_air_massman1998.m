function out = water_vapor_diffusivity_air_massman1998(Tair_K,Tleaf_K,pressure_Pa)
% WATER_VAPOR_DIFFUSIVITY_AIR_MASSMAN1998
%
% Candidate production helper for molecular diffusivity of H2O vapor in air.
%
% SOURCE
% ------
% Massman, W.J. (1998)
% A review of the molecular diffusivities of H2O, CO2, CH4, CO, O3, SO2,
% NH3, N2O, NO, and NO2 in air, O2 and N2 near STP.
% Atmospheric Environment 32(6):1111-1127.
% DOI: 10.1016/S1352-2310(97)00391-9
%
% PARAMETERIZATION
% ----------------
% Dv = D0 * (Tfilm/T0)^b * (P0/P)
%
% D0 = 2.178e-5 m^2 s^-1
% T0 = 273.15 K
% b  = 1.81
% P0 = 101325 Pa
%
% Tfilm = 0.5*(Tair_K + Tleaf_K)

out=struct();
out.success=false;
out.code="UNSET";

out.source="Massman 1998 Atmospheric Environment 32(6):1111-1127";
out.DOI="10.1016/S1352-2310(97)00391-9";

out.referenceDiffusivity_m2_s=2.178e-5;
out.referenceTemperature_K=273.15;
out.referencePressure_Pa=101325;
out.temperatureExponent=1.81;

out.Tair_K=NaN;
out.Tleaf_K=NaN;
out.Tfilm_K=NaN;
out.pressure_Pa=NaN;
out.waterVaporDiffusivity_m2_s=NaN;

vals=[Tair_K,Tleaf_K,pressure_Pa];
if ~(isnumeric(vals) && isreal(vals) && all(isfinite(vals)))
    out.code="INVALID_INPUT";
    return
end

if ~(isscalar(Tair_K) && Tair_K>0)
    out.code="INVALID_TAIR";
    return
end
if ~(isscalar(Tleaf_K) && Tleaf_K>0)
    out.code="INVALID_TLEAF";
    return
end
if ~(isscalar(pressure_Pa) && pressure_Pa>0)
    out.code="INVALID_PRESSURE";
    return
end

Tfilm=0.5*(double(Tair_K)+double(Tleaf_K));

D0=out.referenceDiffusivity_m2_s;
T0=out.referenceTemperature_K;
P0=out.referencePressure_Pa;
b=out.temperatureExponent;

Dv=D0*(Tfilm/T0)^b*(P0/double(pressure_Pa));

if ~(isfinite(Dv) && Dv>0)
    out.code="NONPHYSICAL_RESULT";
    return
end

out.Tair_K=double(Tair_K);
out.Tleaf_K=double(Tleaf_K);
out.Tfilm_K=Tfilm;
out.pressure_Pa=double(pressure_Pa);
out.waterVaporDiffusivity_m2_s=Dv;
out.success=true;
out.code="OK";
end
