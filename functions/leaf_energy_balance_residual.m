function out = leaf_energy_balance_residual( ...
    Tleaf_K,hydraulicE_mol_m2_s,alphaSurface, ...
    shortwaveIncident_W_m2,LWincident_W_m2,epsilonLeaf, ...
    heatArgs)
% LEAF_ENERGY_BALANCE_RESIDUAL
%
% Assemble the steady-state leaf energy-balance residual from
% production components.
%
% Residual:
%   R = Rabs_SW + Rabs_LW - Remit_LW - H_total - lambdaE
%
% INPUTS
% ------
% Tleaf_K
%   scalar leaf temperature [K].
%
% hydraulicE_mol_m2_s
%   authoritative hydro.E [mol H2O m^-2 s^-1].
%
% alphaSurface
%   2-element shortwave absorptance vector:
%       [ABAXIAL, ADAXIAL]
%
% shortwaveIncident_W_m2
%   2-element incident leaf-plane shortwave:
%       [ABAXIAL, ADAXIAL]
%
% LWincident_W_m2
%   2-element incident longwave:
%       [ABAXIAL, ADAXIAL]
%
% epsilonLeaf
%   scalar leaf longwave emissivity.
%
% heatArgs
%   struct containing mixed_surface_heat inputs:
%       uLeaf_m_s
%       forcedCharacteristicDimension_m
%       leafLength_m
%       leafWidth_m
%       inclinationFromHorizontal_deg
%       phi_deg
%       upperAnatomicalSide
%       Tair_K
%       relativeHumidity
%       pressure_Pa
%
%   Air-transport execution modes:
%
%   A) Dynamic Q4dby1 production mode:
%       heatArgs.airTransportMode = "FILM_TEMPERATURE_PRESSURE_DEPENDENT"
%       momentum/heat diffusivities are computed from the current trial Tleaf
%       by air_transport_diffusivities.
%
%   B) Legacy fixed-input mode:
%       momentumDiffusivity_m2_s
%       heatDiffusivity_m2_s
%
% OUTPUT
% ------
% out.success
% out.code
% out.residual_W_m2
% out.RabsSW_W_m2
% out.RabsLW_W_m2
% out.RemitLW_W_m2
% out.Hsurface_W_m2
% out.Htotal_W_m2
% out.lambdaE_W_m2
% out.Tleaf_K
% out.surfaceOrder
% out.mixedHeat
% out.latentHeat


arguments
    Tleaf_K (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
    hydraulicE_mol_m2_s (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative}
    alphaSurface {mustBeNumeric,mustBeReal}
    shortwaveIncident_W_m2 {mustBeNumeric,mustBeReal}
    LWincident_W_m2 {mustBeNumeric,mustBeReal}
    epsilonLeaf (1,1) double {mustBeReal,mustBeFinite}
    heatArgs (1,1) struct
end

surfaceOrder=["ABAXIAL","ADAXIAL"];

alphaSurface=double(alphaSurface(:).');
sw=double(shortwaveIncident_W_m2(:).');
lw=double(LWincident_W_m2(:).');

requireLocal(numel(alphaSurface)==2, ...
    'leaf_energy_balance_residual:AlphaShape', ...
    'alphaSurface must contain [ABAXIAL, ADAXIAL].');

requireLocal(numel(sw)==2, ...
    'leaf_energy_balance_residual:ShortwaveShape', ...
    'shortwaveIncident_W_m2 must contain [ABAXIAL, ADAXIAL].');

requireLocal(numel(lw)==2, ...
    'leaf_energy_balance_residual:LongwaveShape', ...
    'LWincident_W_m2 must contain [ABAXIAL, ADAXIAL].');

requireLocal(all(isfinite(alphaSurface)) && ...
             all(alphaSurface>=0 & alphaSurface<=1), ...
    'leaf_energy_balance_residual:AlphaPhysical', ...
    'alphaSurface must contain finite absorptances in [0,1].');

requireLocal(all(isfinite(sw)) && all(sw>=0), ...
    'leaf_energy_balance_residual:ShortwavePhysical', ...
    'Incident shortwave must be finite and nonnegative.');

requireLocal(all(isfinite(lw)) && all(lw>=0), ...
    'leaf_energy_balance_residual:LongwavePhysical', ...
    'Incident longwave must be finite and nonnegative.');

requireLocal(epsilonLeaf>0 && epsilonLeaf<=1, ...
    'leaf_energy_balance_residual:EpsilonPhysical', ...
    'epsilonLeaf must lie in (0,1].');

requiredHeatFields=[ ...
    "uLeaf_m_s"
    "forcedCharacteristicDimension_m"
    "leafLength_m"
    "leafWidth_m"
    "inclinationFromHorizontal_deg"
    "phi_deg"
    "upperAnatomicalSide"
    "Tair_K"
    "relativeHumidity"
    "pressure_Pa"
    ];

for i=1:numel(requiredHeatFields)
    requireLocal(isfield(heatArgs,requiredHeatFields(i)), ...
        'leaf_energy_balance_residual:MissingHeatArg', ...
        sprintf('heatArgs.%s is required.',requiredHeatFields(i)));
end


airTransportMode="LEGACY_FIXED_INPUT";
if isfield(heatArgs,'airTransportMode')
    airTransportMode=upper(string(heatArgs.airTransportMode));
end

switch airTransportMode
    case "FILM_TEMPERATURE_PRESSURE_DEPENDENT"
        airTransport=air_transport_diffusivities( ...
            heatArgs.Tair_K,Tleaf_K,heatArgs.pressure_Pa);
        momentumDiffusivity_m2_s=airTransport.momentumDiffusivity_m2_s;
        heatDiffusivity_m2_s=airTransport.heatDiffusivity_m2_s;

    case "LEGACY_FIXED_INPUT"
        requireLocal(isfield(heatArgs,'momentumDiffusivity_m2_s'), ...
            'leaf_energy_balance_residual:MissingMomentumDiffusivity', ...
            ['heatArgs.momentumDiffusivity_m2_s is required in legacy ' ...
             'fixed-input mode.']);
        requireLocal(isfield(heatArgs,'heatDiffusivity_m2_s'), ...
            'leaf_energy_balance_residual:MissingHeatDiffusivity', ...
            ['heatArgs.heatDiffusivity_m2_s is required in legacy ' ...
             'fixed-input mode.']);

        momentumDiffusivity_m2_s=double(heatArgs.momentumDiffusivity_m2_s);
        heatDiffusivity_m2_s=double(heatArgs.heatDiffusivity_m2_s);

        requireLocal(isscalar(momentumDiffusivity_m2_s) && ...
            isfinite(momentumDiffusivity_m2_s) && momentumDiffusivity_m2_s>0, ...
            'leaf_energy_balance_residual:InvalidMomentumDiffusivity', ...
            'Legacy momentum diffusivity must be finite, scalar, and positive.');
        requireLocal(isscalar(heatDiffusivity_m2_s) && ...
            isfinite(heatDiffusivity_m2_s) && heatDiffusivity_m2_s>0, ...
            'leaf_energy_balance_residual:InvalidHeatDiffusivity', ...
            'Legacy heat diffusivity must be finite, scalar, and positive.');

        airTransport=struct();
        airTransport.success=true;
        airTransport.code="LEGACY_FIXED_INPUT";
        airTransport.Tfilm_K=0.5*(heatArgs.Tair_K+Tleaf_K);
        airTransport.pressure_kPa=heatArgs.pressure_Pa/1000;
        airTransport.momentumDiffusivity_m2_s=momentumDiffusivity_m2_s;
        airTransport.heatDiffusivity_m2_s=heatDiffusivity_m2_s;

    otherwise
        error('leaf_energy_balance_residual:InvalidAirTransportMode', ...
            'Unsupported heatArgs.airTransportMode: %s',airTransportMode);
end

% Shortwave absorption.
RabsSW_surface=alphaSurface.*sw;
RabsSW=sum(RabsSW_surface);

% Longwave absorption.
RabsLW_surface=epsilonLeaf.*lw;
RabsLW=sum(RabsLW_surface);

% Two-surface longwave emission.
sigma=5.670374419e-8;
Remit_surface=epsilonLeaf.*sigma.*Tleaf_K.^4.*[1 1];
Remit=sum(Remit_surface);

% Sensible heat from existing production helper.
H=mixed_surface_heat( ...
    heatArgs.uLeaf_m_s, ...
    heatArgs.forcedCharacteristicDimension_m, ...
    heatArgs.leafLength_m, ...
    heatArgs.leafWidth_m, ...
    heatArgs.inclinationFromHorizontal_deg, ...
    heatArgs.phi_deg, ...
    heatArgs.upperAnatomicalSide, ...
    momentumDiffusivity_m2_s, ...
    heatDiffusivity_m2_s, ...
    heatArgs.Tair_K, ...
    Tleaf_K, ...
    heatArgs.relativeHumidity, ...
    heatArgs.pressure_Pa);

if ~isfield(H,'success') || ~logical(H.success)
    code="UNKNOWN_MIXED_HEAT_FAILURE";
    if isfield(H,'code')
        code=string(H.code);
    end
    error('leaf_energy_balance_residual:MixedHeatFailure', ...
        'mixed_surface_heat failed: %s',code);
end

requireLocal(isfield(H,'Hsurface_W_m2'), ...
    'leaf_energy_balance_residual:MixedHeatFieldMissing', ...
    'mixed_surface_heat lacks Hsurface_W_m2.');

Hsurface=double(H.Hsurface_W_m2(:).');
requireLocal(numel(Hsurface)==2 && all(isfinite(Hsurface)), ...
    'leaf_energy_balance_residual:MixedHeatShape', ...
    'Hsurface_W_m2 must contain two finite anatomical-surface values.');

Htotal=sum(Hsurface);

% Latent heat from authoritative hydraulic E.
L=leaf_latent_heat_flux(hydraulicE_mol_m2_s,Tleaf_K);
lambdaE=double(L.lambdaE_W_m2);

residual=RabsSW+RabsLW-Remit-Htotal-lambdaE;

out=struct();
out.success=isfinite(residual);
out.code="OK";
out.surfaceOrder=surfaceOrder;
out.Tleaf_K=Tleaf_K;

out.RabsSW_surface_W_m2=RabsSW_surface;
out.RabsSW_W_m2=RabsSW;

out.RabsLW_surface_W_m2=RabsLW_surface;
out.RabsLW_W_m2=RabsLW;

out.RemitLW_surface_W_m2=Remit_surface;
out.RemitLW_W_m2=Remit;

out.Hsurface_W_m2=Hsurface;
out.Htotal_W_m2=Htotal;

out.lambdaE_W_m2=lambdaE;
out.residual_W_m2=residual;

out.mixedHeat=H;
out.latentHeat=L;
out.airTransport=airTransport;
out.airTransportMode=airTransportMode;

out.hydraulicEAuthoritative=true;
out.hydraulicFeedback=false;
out.rootSolvePerformed=false;
end

function requireLocal(condition,id,message)
if ~condition
    error(id,'%s',message);
end
end
