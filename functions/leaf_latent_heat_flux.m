function out = leaf_latent_heat_flux(E_mol_m2_s,Tleaf_K)
% LEAF_LATENT_HEAT_FLUX
%
% Convert molar transpiration flux to latent energy flux.
%
% INPUTS
% ------
% E_mol_m2_s
%   Water-vapour flux on the model one-sided projected lamina-area basis
%   [mol H2O m^-2 s^-1]. Scalar or numeric array. Must be finite and >= 0.
%
% Tleaf_K
%   Leaf temperature [K]. Scalar or numeric array.
%
% OUTPUT
% ------
% out.success
% out.code
% out.E_mol_m2_s
% out.Tleaf_K
% out.Tleaf_C
% out.Lv_J_kg
% out.massFlux_kg_m2_s
% out.lambdaE_W_m2
% out.Mwater_kg_mol
% out.sourceFormula
%
% LATENT-HEAT LAW
% ----------------
% FAO-56 Annex 3 / Harrison (1963):
%
%   Lv_MJ_kg = 2.501 - 2.361e-3*T_C
%
% implemented as:
%
%   Lv_J_kg = (2.501 - 2.361e-3*T_C)*1e6
%
% and:
%
%   lambdaE = E * Mwater * Lv
%
% where:
%   Mwater = 0.01801528 kg mol^-1


arguments
    E_mol_m2_s {mustBeNumeric,mustBeReal}
    Tleaf_K {mustBeNumeric,mustBeReal}
end

E=double(E_mol_m2_s);
Tk=double(Tleaf_K);

if any(~isfinite(E),'all')
    error('leaf_latent_heat_flux:NonfiniteE', ...
        'E_mol_m2_s must be finite.');
end

if any(E<0,'all')
    error('leaf_latent_heat_flux:NegativeE', ...
        'E_mol_m2_s must be nonnegative.');
end

if any(~isfinite(Tk),'all') || any(Tk<=0,'all')
    error('leaf_latent_heat_flux:InvalidTemperature', ...
        'Tleaf_K must be finite and positive.');
end

[E,Tk]=compatibleArraysLocal(E,Tk);

Tc=Tk-273.15;

if any(Tc < -100 | Tc > 100,'all')
    error('leaf_latent_heat_flux:TemperatureOutsideContract', ...
        'Tleaf must lie within [-100,100] degC for this latent-heat law.');
end

Mwater=0.01801528; % kg mol^-1
Lv=(2.501 - 2.361e-3.*Tc).*1e6; % J kg^-1

if any(~isfinite(Lv) | Lv<=0,'all')
    error('leaf_latent_heat_flux:InvalidLatentHeat', ...
        'Computed latent heat must remain finite and positive.');
end

massFlux=E.*Mwater;
lambdaE=massFlux.*Lv;

out=struct();
out.success=true;
out.code="OK";
out.E_mol_m2_s=E;
out.Tleaf_K=Tk;
out.Tleaf_C=Tc;
out.Lv_J_kg=Lv;
out.massFlux_kg_m2_s=massFlux;
out.lambdaE_W_m2=lambdaE;
out.Mwater_kg_mol=Mwater;
out.sourceFormula="Lv_MJ_kg = 2.501 - 2.361e-3*T_C";
out.signConvention="POSITIVE_LAMBDAE_IS_ENERGY_LEAVING_LEAF";
out.areaConvention="ONE_SIDED_PROJECTED_LAMINA_AREA";
out.hydraulicEAuthoritative=true;
end

function [a,b]=compatibleArraysLocal(a,b)
if isscalar(a) && ~isscalar(b)
    a=repmat(a,size(b));
elseif isscalar(b) && ~isscalar(a)
    b=repmat(b,size(a));
elseif ~isequal(size(a),size(b))
    error('leaf_latent_heat_flux:SizeMismatch', ...
        'Inputs must be scalar or have identical sizes.');
end
end
