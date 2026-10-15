function G = graefe2022_free_surface_scalar( ...
    leafLength_m,leafWidth_m,inclinationFromHorizontal_deg,phi_deg, ...
    nu_m2_s,alpha_m2_s,scalarDiffusivity_m2_s, ...
    Tair_K,Tleaf_K,RH,pressure_Pa)
% GRAEFE2022_FREE_SURFACE_SCALAR
%
% Dormant generic helper for free-convective scalar
% boundary-layer conductance over a finite elliptical leaf.
%
% SOURCE BASIS
% ------------
% Graefe et al. (2022), Section 4.5:
%
%   gb = (Ds/kappa)^k * kappa * NuF
%
% with:
%   kc    = 1.00  conduction
%   kfree = 0.75  laminar free convection
%
% and side-specific:
%
%   k = kc - (kc-kfree) * (1 - NuC/max(NuF,NuC))
%
% HEAT-PHYSICS DEPENDENCY
% -----------------------
% The finite-leaf free-convection physics is delegated to:
%
%   graefe2022_free_surface_heat
%
% OUTPUT SIDE ORDER
% -----------------
% [LOWER, UPPER] geometric sides.
%
% PROJECT MOLAR CONVERSION
% ------------------------
%   c_air = pressure_Pa/(R*Tair_K)
%   gb_molar = gb_velocity*c_air
%
% This is the project vapor-gb unit convention.


validateattributes(scalarDiffusivity_m2_s,{'numeric'}, ...
    {'scalar','real','finite','positive'});

H = graefe2022_free_surface_heat( ...
    leafLength_m,leafWidth_m,inclinationFromHorizontal_deg,phi_deg, ...
    nu_m2_s,alpha_m2_s,Tair_K,Tleaf_K,RH,pressure_Pa);

assert(H.success, ...
    'graefe2022_free_surface_scalar:HeatKernelFailure', ...
    'graefe2022_free_surface_heat did not return success.');

NuF = H.transport.NuF_m1;
NuC = H.transport.NuConduction_m1;

conductionFraction = NuC ./ max(NuF,NuC);
conductionFraction = min(max(conductionFraction,0),1);

kConduction = 1.0;
kLaminarFree = 0.75;
kExponent = kConduction - (kConduction-kLaminarFree) .* ...
    (1-conductionFraction);

diffusivityRatio = scalarDiffusivity_m2_s/alpha_m2_s;

gbVelocity_m_s = diffusivityRatio.^kExponent .* ...
    alpha_m2_s .* NuF;

Rgas_J_mol_K = 8.3144598;
airMolarDensity_mol_m3 = pressure_Pa/(Rgas_J_mol_K*Tair_K);
gbMolar_mol_m2_s = gbVelocity_m_s .* airMolarDensity_mol_m3;

G = struct();
G.success = true;
G.code = "OK";
G.surfaceOrder = ["LOWER","UPPER"];

G.scalar = struct();
G.scalar.diffusivity_m2_s = scalarDiffusivity_m2_s;
G.scalar.thermalDiffusivity_m2_s = alpha_m2_s;
G.scalar.diffusivityRatio = diffusivityRatio;
G.scalar.kConduction = kConduction;
G.scalar.kLaminarFree = kLaminarFree;
G.scalar.kExponent = kExponent;
G.scalar.conductionFraction = conductionFraction;

G.transport = struct();
G.transport.NuF_m1 = NuF;
G.transport.NuConduction_m1 = NuC;
G.transport.gbVelocity_m_s = gbVelocity_m_s;
G.transport.airMolarDensity_mol_m3 = airMolarDensity_mol_m3;
G.transport.gbMolar_mol_m2_s = gbMolar_mol_m2_s;

G.heatKernel = H;

G.provenance = struct();
G.provenance.source = "Graefe_et_al_2022_Section_4_5";
G.provenance.heatKernel = "graefe2022_free_surface_heat";
G.provenance.productionPromotionStage = "Q3f3j105";
G.provenance.productionImplementationStage = "Q3f3j106";
G.provenance.scope = "FREE_SCALAR_TRANSFER_ONLY";
G.provenance.molarConversion = "PROJECT_P_OVER_R_TAIR";
G.provenance.B2020ConditioningEmbedded = false;
G.provenance.anatomicalMappingEmbedded = false;
end
