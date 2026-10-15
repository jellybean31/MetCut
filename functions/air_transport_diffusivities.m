function out = air_transport_diffusivities(Tair_K,Tleaf_K,pressure_Pa)
% AIR_TRANSPORT_DIFFUSIVITIES
%
% Production helper for state-dependent air momentum and heat
% diffusivities.
%
% Tfilm_K = 0.5*(Tair_K + Tleaf_K)
% D = D0 * (Tfilm_K/273.15)^1.75 * (101.3246/pressure_kPa)
%
% Reference coefficients:
%   momentum diffusivity Dm0 = 13.3e-6 m^2 s^-1
%   heat diffusivity     Dh0 = 19.0e-6 m^2 s^-1

arguments
    Tair_K (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
    Tleaf_K (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
    pressure_Pa (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
end

T0_K = 273.15;
P0_kPa = 101.3246;
eT = 1.75;
Dm0 = 13.3e-6;
Dh0 = 19.0e-6;

Tfilm_K = 0.5*(Tair_K + Tleaf_K);
pressure_kPa = pressure_Pa/1000;

scale = (Tfilm_K/T0_K)^eT * (P0_kPa/pressure_kPa);
Dm = Dm0*scale;
Dh = Dh0*scale;

if ~(isfinite(Dm) && Dm>0 && isfinite(Dh) && Dh>0)
    error('air_transport_diffusivities:InvalidOutput', ...
        'Computed diffusivities must be finite and positive.');
end

out = struct();
out.success = true;
out.code = "OK";
out.Tfilm_K = Tfilm_K;
out.pressure_kPa = pressure_kPa;
out.momentumDiffusivity_m2_s = Dm;
out.heatDiffusivity_m2_s = Dh;
out.referenceTemperature_K = T0_K;
out.referencePressure_kPa = P0_kPa;
out.temperatureExponent = eT;
out.momentumReference_m2_s = Dm0;
out.heatReference_m2_s = Dh0;
out.scope = "AIR_TRANSPORT_PROPERTIES_ONLY";
out.hydraulicFeedback = false;
end
