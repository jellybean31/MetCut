function out = mixed_surface_scalar( ...
    uLeaf_m_s,forcedCharacteristicDimension_m, ...
    leafLength_m,leafWidth_m,inclinationFromHorizontal_deg,phi_deg, ...
    upperAnatomicalSide,momentumDiffusivity_m2_s,heatDiffusivity_m2_s, ...
    waterVaporDiffusivity_m2_s,Tair_K,Tleaf_K,relativeHumidity,pressure_Pa)
% MIXED_SURFACE_SCALAR
%
% CANDIDATE generic mixed forced/free vapor boundary-layer conductance helper.
%
% COMPONENTS
% ----------
% Forced scalar transfer:
%   forced_surface_gb.m
%   native/output order [ABAXIAL, ADAXIAL]
%
% Free scalar transfer:
%   graefe2022_free_surface_scalar.m
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
% OUTPUT ORDER
% ------------
% [ABAXIAL, ADAXIAL]


out=struct();
out.success=false;
out.code="UNSET";
out.scope="VAPOR_SCALAR_MIXED_CONVECTION_CANDIDATE";
out.surfaceOrder=["ABAXIAL","ADAXIAL"];
out.mixingBasis="VELOCITY_CONDUCTANCE_M_S";
out.mixingExponent=3.5;
out.upperAnatomicalSide=string(upperAnatomicalSide);
out.inclinationFromHorizontal_deg=inclinationFromHorizontal_deg;

out.forced=struct();
out.freeGeometric=struct();
out.freeMappingVelocity=struct();
out.freeMappingMolar=struct();

out.gbForcedVelocity_m_s=[NaN NaN];
out.gbFreeVelocity_m_s=[NaN NaN];
out.gbMixedVelocity_m_s=[NaN NaN];

out.gbForced_mol_m2_s=[NaN NaN];
out.gbFree_mol_m2_s=[NaN NaN];
out.gbMixed_mol_m2_s=[NaN NaN];

out.airMolarDensityForced_mol_m3=NaN;
out.airMolarDensityFree_mol_m3=NaN;
out.airMolarDensity_mol_m3=NaN;
out.airMolarDensityResidual_mol_m3=NaN;

% Basic scalar validations kept local so component failures remain explicit.
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

if ~(isnumeric(waterVaporDiffusivity_m2_s) && ...
        isreal(waterVaporDiffusivity_m2_s) && ...
        isscalar(waterVaporDiffusivity_m2_s) && ...
        isfinite(waterVaporDiffusivity_m2_s) && ...
        waterVaporDiffusivity_m2_s>0)
    out.code="INVALID_WATER_VAPOR_DIFFUSIVITY";
    return
end

F=forced_surface_gb( ...
    uLeaf_m_s,forcedCharacteristicDimension_m, ...
    momentumDiffusivity_m2_s,heatDiffusivity_m2_s, ...
    waterVaporDiffusivity_m2_s,Tair_K,pressure_Pa);

out.forced=F;

if ~isfield(F,'success') || ~F.success
    out.code="FORCED_COMPONENT_FAILED";
    if isfield(F,'code')
        out.componentCode=string(F.code);
    end
    return
end

G=graefe2022_free_surface_scalar( ...
    leafLength_m,leafWidth_m,inclinationFromHorizontal_deg,phi_deg, ...
    momentumDiffusivity_m2_s,heatDiffusivity_m2_s, ...
    waterVaporDiffusivity_m2_s, ...
    Tair_K,Tleaf_K,relativeHumidity,pressure_Pa);

out.freeGeometric=G;

if ~isfield(G,'success') || ~G.success
    out.code="FREE_COMPONENT_FAILED";
    if isfield(G,'code')
        out.componentCode=string(G.code);
    end
    return
end

if ~(isfield(F,'gbVelocity_m_s') && ...
        isfield(F,'gb_mol_m2_s') && ...
        isfield(F,'airMolarDensity_mol_m3'))
    out.code="FORCED_OUTPUT_SCHEMA_MISMATCH";
    return
end

if ~(isfield(G,'transport') && ...
        isfield(G.transport,'gbVelocity_m_s') && ...
        isfield(G.transport,'gbMolar_mol_m2_s') && ...
        isfield(G.transport,'airMolarDensity_mol_m3'))
    out.code="FREE_OUTPUT_SCHEMA_MISMATCH";
    return
end

gForcedV=reshape(double(F.gbVelocity_m_s),1,2);
gForcedM=reshape(double(F.gb_mol_m2_s),1,2);

gFreeGeomV=reshape(double(G.transport.gbVelocity_m_s),1,2);
gFreeGeomM=reshape(double(G.transport.gbMolar_mol_m2_s),1,2);

MV=map_geometric_two_surface_to_anatomical( ...
    gFreeGeomV,upperAnatomicalSide,inclinationFromHorizontal_deg);
MM=map_geometric_two_surface_to_anatomical( ...
    gFreeGeomM,upperAnatomicalSide,inclinationFromHorizontal_deg);

out.freeMappingVelocity=MV;
out.freeMappingMolar=MM;

if ~MV.success || ~MM.success
    out.code="FREE_SURFACE_MAPPING_FAILED";
    if ~MV.success
        out.componentCode=MV.code;
    else
        out.componentCode=MM.code;
    end
    return
end

gFreeV=reshape(double(MV.anatomicalValues),1,2);
gFreeM=reshape(double(MM.anatomicalValues),1,2);

if any(~isfinite([gForcedV gFreeV gForcedM gFreeM])) || ...
        any([gForcedV gFreeV gForcedM gFreeM] < 0)
    out.code="INVALID_COMPONENT_CONDUCTANCE";
    return
end

cForced=double(F.airMolarDensity_mol_m3);
cFree=double(G.transport.airMolarDensity_mol_m3);

if ~(isscalar(cForced) && isfinite(cForced) && cForced>0 && ...
        isscalar(cFree) && isfinite(cFree) && cFree>0)
    out.code="INVALID_COMPONENT_AIR_MOLAR_DENSITY";
    return
end

densityResidual=abs(cForced-cFree);
densityTol=1e-12*max([abs(cForced),abs(cFree),1]);

if densityResidual>densityTol
    out.code="AIR_MOLAR_DENSITY_CONVENTION_MISMATCH";
    out.airMolarDensityForced_mol_m3=cForced;
    out.airMolarDensityFree_mol_m3=cFree;
    out.airMolarDensityResidual_mol_m3=densityResidual;
    return
end

p=out.mixingExponent;
gMixedV=(gForcedV.^p + gFreeV.^p).^(1/p);

% Common scalar-transfer molar conversion demonstrated in Q4dby23a-v2.
cAir=0.5*(cForced+cFree);
gMixedM=gMixedV*cAir;

% Internal consistency against component-native molar values.
forcedResidual=max(abs(gForcedM-gForcedV*cAir));
freeResidual=max(abs(gFreeM-gFreeV*cAir));
conversionTol=1e-10;

if forcedResidual>conversionTol || freeResidual>conversionTol
    out.code="COMPONENT_MOLAR_CONVERSION_MISMATCH";
    out.forcedMolarConversionResidual=forcedResidual;
    out.freeMolarConversionResidual=freeResidual;
    return
end

out.gbForcedVelocity_m_s=gForcedV;
out.gbFreeVelocity_m_s=gFreeV;
out.gbMixedVelocity_m_s=gMixedV;

out.gbForced_mol_m2_s=gForcedM;
out.gbFree_mol_m2_s=gFreeM;
out.gbMixed_mol_m2_s=gMixedM;

out.airMolarDensityForced_mol_m3=cForced;
out.airMolarDensityFree_mol_m3=cFree;
out.airMolarDensity_mol_m3=cAir;
out.airMolarDensityResidual_mol_m3=densityResidual;

out.forcedMolarConversionResidual=forcedResidual;
out.freeMolarConversionResidual=freeResidual;

out.provenance=struct();
out.provenance.implementationRole="CANDIDATE_GENERIC";
out.provenance.forcedHelper="forced_surface_gb";
out.provenance.freeHelper="graefe2022_free_surface_scalar";
out.provenance.surfaceMapper="map_geometric_two_surface_to_anatomical";
out.provenance.mixingArchitecture="MATCHES_MIXED_SURFACE_HEAT_VELOCITY_P3P5";
out.provenance.scalarOrVaporExtension=true;
out.provenance.heatConductanceReusedAsVapor=false;
out.provenance.B2020DefaultFacing=false;
out.provenance.hydraulicFeedback=false;
out.provenance.numericalTleafSolver=false;
out.provenance.developmentStage="Q4dby23b";

out.success=true;
out.code="OK";
end
