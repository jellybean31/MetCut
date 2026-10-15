
function T = wax_accumulation_default_traits()
%WAX_ACCUMULATION_DEFAULT_TRAITS
% Default/development configuration for standing leaf-wax accumulation.
%
% The accumulation layer is intentionally separate from instantaneous wax
% biosynthesis. The existing biochemical output remains the isotope
% composition of newly synthesized wax.
%
% Mass units are arbitrary but internally consistent. By default, peak
% juvenile production is 1 mass unit d^-1 for 14 d, so the reference mature
% wax mass is 14 units.
%
% PRODUCTION MODES
%   continuous
%   developmentGated
%   developmentPlusMaintenance
%
% developmentPlusMaintenance uses a mature replacement fraction per day.
% Mature production is replacementFractionPerDay * referenceMatureMass.
% If lossRatePerDay equals replacementFractionPerDay, mature wax mass is
% approximately stationary after leaf development.

T=struct;

T.productionMode="developmentPlusMaintenance";

T.developmentWindowDays=14;
T.peakProductionRate_perDay=1;

% Developmental shape:
%   flat  -> constant juvenile production
%   sine  -> smooth zero-at-edges pulse over development window
T.developmentShape="flat";

% Initial standing pool.
T.initialMass=0;
T.initialDeltaD_permil=NaN;

% Mature secondary synthesis / renewal.
% Development value only; empirical fitting/validation should replace it.
T.matureReplacementFractionPerDay=0.005;

% Standing-pool loss/removal. Keeping this equal to the mature replacement
% fraction gives approximately constant mature mass under maintenance mode.
T.lossRatePerDay=0.005;

% Homologue-specific multipliers may be scalar or nHomologue x 1.
T.productionMultiplier=1;
T.lossMultiplier=1;

% Numerical tolerance for defining a non-empty pool.
T.massTolerance=1e-14;

% Standard isotope ratio used only for exact deltaD <-> D atom fraction
% conversion inside this layer.
T.Rstd_VSMOW=155.76e-6;

end
