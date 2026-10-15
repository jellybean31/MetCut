function out = mixed_surface_heat( ...
    uLeaf_m_s,forcedCharacteristicDimension_m, ...
    leafLength_m,leafWidth_m,inclinationFromHorizontal_deg,phi_deg, ...
    upperAnatomicalSide,momentumDiffusivity_m2_s,heatDiffusivity_m2_s, ...
    Tair_K,Tleaf_K,relativeHumidity,pressure_Pa)
% MIXED_SURFACE_HEAT
%
% Generic production sensible-heat mixed-convection helper.
%
% COMPONENTS
% ----------
% Forced convection:
%   forced_surface_heat.m
%   native/output order [ABAXIAL, ADAXIAL]
%
% Free convection:
%   graefe2022_free_surface_heat.m
%   native geometric order [LOWER, UPPER]
%
% Surface mapping:
%   map_geometric_two_surface_to_anatomical.m
%
% MIXING BASIS
% ------------
% Velocity conductance [m s^-1]:
%
%   gMixed = (gForced^p + gFree^p)^(1/p)
%
% with:
%   p = 3.5
%
%
% INPUTS
% ------
% uLeaf_m_s
%   2-element forced-convection leaf-level wind [m s^-1], anatomical order
%   [ABAXIAL, ADAXIAL].
%
% forcedCharacteristicDimension_m
%   2-element forced-convection characteristic dimension [m], anatomical
%   order [ABAXIAL, ADAXIAL].
%
% leafLength_m, leafWidth_m
%   Positive leaf dimensions [m] used by Graefe free convection.
%   leafWidth_m <= leafLength_m.
%
% inclinationFromHorizontal_deg
%   0 = horizontal, 90 = vertical.
%
% phi_deg
%   Graefe streamline-orientation angle, 0..90 deg.
%
% upperAnatomicalSide
%   "ADAXIAL", "ABAXIAL", or "UNRESOLVED".
%   UNRESOLVED is admissible only through the production mapper's exact-
%   vertical geometric-degeneracy rule.
%
% momentumDiffusivity_m2_s
% heatDiffusivity_m2_s
% Tair_K
% Tleaf_K
% relativeHumidity
% pressure_Pa
%   Passed to the production component helpers.
%
% OUTPUT ORDER
% ------------
% [ABAXIAL, ADAXIAL]
%
% H > 0 means sensible heat leaves the leaf.

out=struct();
out.success=false;
out.code="UNSET";
out.scope="SENSIBLE_HEAT_MIXED_CONVECTION_ONLY";
out.surfaceOrder=["ABAXIAL","ADAXIAL"];
out.mixingBasis="VELOCITY_CONDUCTANCE_M_S";
out.mixingExponent=3.5;
out.upperAnatomicalSide=string(upperAnatomicalSide);
out.inclinationFromHorizontal_deg=inclinationFromHorizontal_deg;

out.forced=struct();
out.freeGeometric=struct();
out.freeMapping=struct();

out.gHForced_m_s=[NaN NaN];
out.gHFree_m_s=[NaN NaN];
out.gHMixed_m_s=[NaN NaN];
out.airMolarDensity_mol_m3=NaN;
out.gHMixed_mol_m2_s=[NaN NaN];
out.cpAirMolar_J_mol_K=NaN;
out.Hsurface_W_m2=[NaN NaN];

if ~(isnumeric(leafLength_m) && isreal(leafLength_m) && ...
        isscalar(leafLength_m) && isfinite(leafLength_m) && leafLength_m>0)
    out.code="INVALID_LEAF_LENGTH";
    return
end

if ~(isnumeric(leafWidth_m) && isreal(leafWidth_m) && ...
        isscalar(leafWidth_m) && isfinite(leafWidth_m) && ...
        leafWidth_m>0 && leafWidth_m<=leafLength_m)
    out.code="INVALID_LEAF_WIDTH";
    return
end

if ~(isnumeric(phi_deg) && isreal(phi_deg) && isscalar(phi_deg) && ...
        isfinite(phi_deg) && phi_deg>=0 && phi_deg<=90)
    out.code="INVALID_PHI";
    return
end

if ~(isnumeric(relativeHumidity) && isreal(relativeHumidity) && ...
        isscalar(relativeHumidity) && isfinite(relativeHumidity) && ...
        relativeHumidity>=0 && relativeHumidity<=1)
    out.code="INVALID_RH";
    return
end

F=forced_surface_heat( ...
    uLeaf_m_s,forcedCharacteristicDimension_m, ...
    momentumDiffusivity_m2_s,heatDiffusivity_m2_s, ...
    Tair_K,Tleaf_K,pressure_Pa);

out.forced=F;

if ~isfield(F,'success') || ~F.success
    out.code="FORCED_COMPONENT_FAILED";
    if isfield(F,'code')
        out.componentCode=string(F.code);
    end
    return
end

G=graefe2022_free_surface_heat( ...
    leafLength_m,leafWidth_m,inclinationFromHorizontal_deg,phi_deg, ...
    momentumDiffusivity_m2_s,heatDiffusivity_m2_s, ...
    Tair_K,Tleaf_K,relativeHumidity,pressure_Pa);

out.freeGeometric=G;

if ~isfield(G,'success') || ~G.success
    out.code="FREE_COMPONENT_FAILED";
    if isfield(G,'code')
        out.componentCode=string(G.code);
    end
    return
end

[gFreeGeometric,freeRecoverCode]=recoverFreeVelocityLocal(G);
if isempty(gFreeGeometric)
    out.code=freeRecoverCode;
    return
end

M=map_geometric_two_surface_to_anatomical( ...
    gFreeGeometric,upperAnatomicalSide,inclinationFromHorizontal_deg);

out.freeMapping=M;

if ~M.success
    out.code="FREE_SURFACE_MAPPING_FAILED";
    out.componentCode=M.code;
    return
end

gForced=reshape(double(F.gHVelocity_m_s),1,2);
gFree=reshape(double(M.anatomicalValues),1,2);

if any(gForced<0) || any(gFree<0)
    out.code="NEGATIVE_COMPONENT_CONDUCTANCE";
    return
end

p=out.mixingExponent;
gMixed=(gForced.^p + gFree.^p).^(1/p);

if ~isfield(F,'airMolarDensity_mol_m3') || ...
        ~isfield(F,'cpAirMolar_J_mol_K')
    out.code="FORCED_THERMODYNAMIC_FIELDS_MISSING";
    return
end

cAir=double(F.airMolarDensity_mol_m3);
cp=double(F.cpAirMolar_J_mol_K);

gMolar=gMixed*cAir;
H=cp*gMolar*(Tleaf_K-Tair_K);

out.gHForced_m_s=gForced;
out.gHFree_m_s=gFree;
out.gHMixed_m_s=gMixed;
out.airMolarDensity_mol_m3=cAir;
out.gHMixed_mol_m2_s=gMolar;
out.cpAirMolar_J_mol_K=cp;
out.Hsurface_W_m2=H;

out.provenance=struct();
out.provenance.implementationRole="PRODUCTION_GENERIC";
out.provenance.forcedHelper="forced_surface_heat";
out.provenance.freeHelper="graefe2022_free_surface_heat";
out.provenance.surfaceMapper="map_geometric_two_surface_to_anatomical";
out.provenance.rawNuCombinationUsed=false;
out.provenance.scalarOrVaporExtension=false;
out.provenance.B2020DefaultFacing=false;
out.provenance.hydraulicFeedback=false;
out.provenance.numericalTleafSolver=false;

out.success=true;
out.code="OK";
end

function [g,code] = recoverFreeVelocityLocal(G)
g=[];
code="FREE_VELOCITY_FIELD_NOT_FOUND";

if isfield(G,'transport')
    T=G.transport;
    names={'gHVelocity_m_s','gH_m_s','conductanceVelocity_m_s'};
    for i=1:numel(names)
        if isfield(T,names{i})
            g=double(T.(names{i}));
            break
        end
    end
end

if isempty(g) && isfield(G,'gHVelocity_m_s')
    g=double(G.gHVelocity_m_s);
end

if isempty(g)
    return
end

if ~(isnumeric(g) && isreal(g) && numel(g)==2 && all(isfinite(g(:))))
    g=[];
    code="INVALID_FREE_VELOCITY_FIELD";
    return
end

g=reshape(g,1,2);
code="OK";
end
