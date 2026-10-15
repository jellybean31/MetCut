function G = graefe2022_free_surface_heat( ...
    leafLength_m,leafWidth_m,inclinationFromHorizontal_deg,phi_deg, ...
    nu_m2_s,alpha_m2_s,Tair_K,Tleaf_K,RH,pressure_Pa)
% GRAEFE2022_FREE_SURFACE_HEAT
%
% Dormant generic production helper for free sensible-heat convection from a
% finite elliptical leaf using the Graefe et al. (2022) formulation, with the
% project virtual-temperature buoyancy.
%
% OUTPUT SIDE ORDER
% -----------------
% [LOWER, UPPER] geometric sides.
%
%
% ORIENTATION
% -----------
% inclinationFromHorizontal_deg:
%   0 = horizontal
%   90 = vertical
%
% phi_deg:
%   streamline direction relative to the ellipse minor axis over the unique
%   symmetry domain [0,90] deg.
%
% BUOYANCY / HEAT-FLUX SEMANTICS
% ------------------------------
% |DeltaTv|:
%   free-convection magnitude
%
% sign(DeltaTv):
%   heated-like / cooled-like geometric branch
%
% sign(actual DeltaT):
%   sensible-heat flux direction; returned but not multiplied into gH here.
%
% REFERENCE AREA / GEOMETRY
% -------------------------
% leafLength_m and leafWidth_m are full ellipse major/minor dimensions.
% One-sided area S = pi*L*W/4.


validateattributes(leafLength_m,{'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(leafWidth_m,{'numeric'}, ...
    {'scalar','real','finite','positive'});
assert(leafLength_m>=leafWidth_m, ...
    'graefe2022_free_surface_heat:LengthWidthOrder', ...
    'leafLength_m must be >= leafWidth_m.');

validateattributes(inclinationFromHorizontal_deg,{'numeric'}, ...
    {'scalar','real','finite','>=',0,'<=',90});
validateattributes(phi_deg,{'numeric'}, ...
    {'scalar','real','finite','>=',0,'<=',90});
validateattributes(nu_m2_s,{'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(alpha_m2_s,{'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(Tair_K,{'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(Tleaf_K,{'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(RH,{'numeric'}, ...
    {'scalar','real','finite','>=',0,'<=',1});
validateattributes(pressure_Pa,{'numeric'}, ...
    {'scalar','real','finite','positive'});

L = leafLength_m;
W = leafWidth_m;
AR = L/W;
S = 0.25*pi*L*W;

% Graefe Eq. (8): ellipse dissection-index approximation.
DI = 1.112*((L+W)/sqrt(L*W))^0.783 - 0.913;

% Graefe Eq. (7), written through DI.
Lh = 0.5*sqrt(S)/(DI*sqrt(pi));

% Graefe Eq. (9): inclined characteristic length for n=3*nRa=0.75.
nInclined = 0.75;
phi = deg2rad(phi_deg);
u = W/2;
v = L/2;

% Appendix-A directional mean-chord expression recovered/audited in j89b.
LCo = pi*u*v/(2*sqrt((u*sin(phi))^2 + (v*cos(phi))^2));
Li = (0.04511*nInclined + 1.006)*LCo;

% Exact project virtual-temperature buoyancy treatment.
epsilonWater = 0.622;
eAir_Pa = svpPaLocal(Tair_K-273.15)*RH;
eLeaf_Pa = svpPaLocal(Tleaf_K-273.15);

assert(eAir_Pa>=0 && eAir_Pa<pressure_Pa, ...
    'graefe2022_free_surface_heat:AirVaporPressure', ...
    'Ambient vapor pressure must satisfy 0 <= eAir < pressure.');
assert(eLeaf_Pa>0 && eLeaf_Pa<pressure_Pa, ...
    'graefe2022_free_surface_heat:LeafVaporPressure', ...
    'Saturated leaf vapor pressure must satisfy 0 < eLeaf < pressure.');

denAir = 1-(1-epsilonWater)*eAir_Pa/pressure_Pa;
denLeaf = 1-(1-epsilonWater)*eLeaf_Pa/pressure_Pa;
assert(denAir>0 && denLeaf>0, ...
    'graefe2022_free_surface_heat:VirtualTemperatureDenominator', ...
    'Virtual-temperature denominator must remain positive.');

TvAir_K = Tair_K/denAir;
TvLeaf_K = Tleaf_K/denLeaf;
deltaTv_K = TvLeaf_K-TvAir_K;
deltaT_K = Tleaf_K-Tair_K;

Pr = nu_m2_s/alpha_m2_s;
g = 9.80665;

% Project-adapted Graefe buoyancy coefficient; dimensions m^-3.
RaTv_m3inv = g*abs(deltaTv_K)*Pr/(Tair_K*nu_m2_s^2);

% Graefe fitted parameters, Table 3.
p1 = 6;
p2 = 1.42;
p3 = 0.513;
p4 = 0.4;
p5 = 0.202;
p6 = 0.492;

s = sind(inclinationFromHorizontal_deg);
c = cosd(inclinationFromHorizontal_deg);

% Heated-like RHS components.
NuInclinedUpperHeated = p3*Li^(-0.25)*s^p4*RaTv_m3inv^0.25;
NuInclinedLowerHeated = p3*Li^(-0.25)*s^0.25*RaTv_m3inv^0.25;

NuHorizontalUpperHeated = p5*Lh^(-0.25)*c^0.25*RaTv_m3inv^0.25;
NuHorizontalLowerHeated = p6*Lh^(-0.4)*c^0.20*RaTv_m3inv^0.20;

% Graefe conduction/diffusion lower limit.
NuConduction_m1 = 5.3*sin(0.0338*AR + 0.607)/sqrt(2*S);

tolTv = 64*eps(max([abs(TvAir_K),abs(TvLeaf_K),1]));

if abs(deltaTv_K)<=tolTv
    direction = "NEUTRAL";
    NuInclinedUpper = 0;
    NuInclinedLower = 0;
    NuHorizontalUpper = 0;
    NuHorizontalLower = 0;
elseif deltaTv_K>0
    direction = "HEATED_LIKE";
    NuInclinedUpper = NuInclinedUpperHeated;
    NuInclinedLower = NuInclinedLowerHeated;
    NuHorizontalUpper = NuHorizontalUpperHeated;
    NuHorizontalLower = NuHorizontalLowerHeated;
else
    direction = "COOLED_LIKE";
    % Graefe cooled-leaf rule: exchange upper/lower RHS.
    NuInclinedUpper = NuInclinedLowerHeated;
    NuInclinedLower = NuInclinedUpperHeated;
    NuHorizontalUpper = NuHorizontalLowerHeated;
    NuHorizontalLower = NuHorizontalUpperHeated;
end

NuThinUpper = (NuInclinedUpper^p1 + NuHorizontalUpper^p1)^(1/p1);
NuThinLower = (NuInclinedLower^p1 + NuHorizontalLower^p1)^(1/p1);

NuFUpper_m1 = (NuThinUpper^p2 + NuConduction_m1^p2)^(1/p2);
NuFLower_m1 = (NuThinLower^p2 + NuConduction_m1^p2)^(1/p2);

gHUpper_m_s = alpha_m2_s*NuFUpper_m1;
gHLower_m_s = alpha_m2_s*NuFLower_m1;

G = struct();
G.success = true;
G.code = "OK";
G.surfaceOrder = ["LOWER","UPPER"];

G.geometry = struct();
G.geometry.leafLength_m = L;
G.geometry.leafWidth_m = W;
G.geometry.areaOneSided_m2 = S;
G.geometry.aspectRatio = AR;
G.geometry.dissectionIndex = DI;
G.geometry.Lh_m = Lh;
G.geometry.LCo_m = LCo;
G.geometry.Li_m = Li;
G.geometry.inclinationFromHorizontal_deg = inclinationFromHorizontal_deg;
G.geometry.phi_deg = phi_deg;

G.thermo = struct();
G.thermo.Tair_K = Tair_K;
G.thermo.Tleaf_K = Tleaf_K;
G.thermo.deltaT_K = deltaT_K;
G.thermo.RH = RH;
G.thermo.pressure_Pa = pressure_Pa;
G.thermo.vaporPressureAir_Pa = eAir_Pa;
G.thermo.vaporPressureLeaf_Pa = eLeaf_Pa;
G.thermo.TvAir_K = TvAir_K;
G.thermo.TvLeaf_K = TvLeaf_K;
G.thermo.deltaTv_K = deltaTv_K;
G.thermo.freeConvectionDirection = direction;
G.thermo.sensibleHeatFluxSign = sign(deltaT_K);

G.transport = struct();
G.transport.Pr = Pr;
G.transport.RaTv_m3inv = RaTv_m3inv;
G.transport.NuInclined_m1 = [NuInclinedLower NuInclinedUpper];
G.transport.NuHorizontal_m1 = [NuHorizontalLower NuHorizontalUpper];
G.transport.NuThin_m1 = [NuThinLower NuThinUpper];
G.transport.NuConduction_m1 = NuConduction_m1;
G.transport.NuF_m1 = [NuFLower_m1 NuFUpper_m1];
G.transport.gHVelocity_m_s = [gHLower_m_s gHUpper_m_s];
G.transport.conductionFraction = ...
    [NuConduction_m1/NuFLower_m1 NuConduction_m1/NuFUpper_m1];

G.provenance = struct();
G.provenance.source = "Graefe_et_al_2022";
G.provenance.productionPromotionStage = "Q3f3j98";
G.provenance.productionImplementationStage = "Q3f3j99";
G.provenance.scope = "FREE_SENSIBLE_HEAT_ONLY";
G.provenance.orientationConditioning = "EXPLICIT_PHI_INPUT";
G.provenance.B2020ConditioningEmbedded = false;
end

function esPa = svpPaLocal(TC)
esPa = 1000*0.61078.*exp((17.2694.*TC)./(TC+237.29));
end
