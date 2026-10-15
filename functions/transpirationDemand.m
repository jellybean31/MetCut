function [E, gs] = transpirationDemand( ...
    psiLeaf, ...
    met, ...
    traits, ...
    t)
%TRANSPIRATIONDEMAND
%
% Stomatal conductance and transpiration demand.
%
% Supports scalar or time-dependent meteorological forcing.
%
% INPUT
%
%   psiLeaf
%       leaf water potential, MPa
%
%   met
%       meteorological forcing structure
%
%   traits
%       plant hydraulic / stomatal traits
%
%   t
%       optional timestep index when met fields are vectors
%
%
% Stomatal model:
%
%   gs =
%       gsMax
%       * fPsi
%       * fLight
%
% where:
%
%   fPsi
%       hydraulic stomatal response
%
%   fLight
%       environmental/light multiplier
%
%
% Backward compatibility:
%
%   If met.lightFraction is absent:
%
%       fLight = 1
%
% and the previous hydraulic stomatal formulation is recovered.


%% ========================================================================
% Optional timestep
% ========================================================================

if nargin < 4

    t = ...
        [];

end


%% ========================================================================
% Meteorological forcing
% ========================================================================

Tleaf = ...
    getTimeValue( ...
        met.TleafC, ...
        t, ...
        'TleafC');


Tair = ...
    getTimeValue( ...
        met.TairC, ...
        t, ...
        'TairC');


RH = ...
    getTimeValue( ...
        met.RH, ...
        t, ...
        'RH');


if isfield(met,'pressureKPa')

    pressure = ...
        getTimeValue( ...
            met.pressureKPa, ...
            t, ...
            'pressureKPa');

else

    pressure = ...
        101.325;

end


if isfield(met,'gb_mol_m2_s')

    gb = ...
        getTimeValue( ...
            met.gb_mol_m2_s, ...
            t, ...
            'gb_mol_m2_s');

else

    gb = ...
        1.0;

end


%% ========================================================================
% Hydraulic stomatal response
% ========================================================================

x = ...
    abs(psiLeaf);


x50 = ...
    abs(traits.psi50Stomata);


fPsi = ...
    1 ...
    ./ ( ...
        1 ...
        + exp( ...
            traits.stomatalSlope ...
            .* (x - x50)));


%% ========================================================================
% Environmental / light response
%
% If no light forcing is supplied, preserve historical behavior.
% ========================================================================

if isfield(met,'lightFraction')

    fLight = ...
        getTimeValue( ...
            met.lightFraction, ...
            t, ...
            'lightFraction');


    fLight = ...
        min( ...
            max(fLight,0), ...
            1);

else

    fLight = ...
        1;

end


%% ========================================================================
% Stomatal conductance
% ========================================================================

gs = ...
    traits.gsMax ...
    .* fPsi ...
    .* fLight;


%% ========================================================================
% Vapor-pressure gradient
% ========================================================================

ei = ...
    saturationVaporPressureKPa( ...
        Tleaf);


ea = ...
    RH ...
    .* saturationVaporPressureKPa( ...
        Tair);


%% ========================================================================
% Total vapor conductance
% ========================================================================

% Handle full stomatal closure explicitly so that gs = 0 does not
% generate numerical warnings through 1/gs.

if gs <= 0

    gw = ...
        0;

else

    gw = ...
        1 ...
        ./ ( ...
            1 ./ gs ...
            + 1 ./ gb);

end


%% ========================================================================
% Transpiration
% ========================================================================

E = ...
    gw ...
    .* max( ...
        ei - ea, ...
        0) ...
    ./ pressure;

end


%% =========================================================================
function x = getTimeValue(x,t,name)

if isscalar(x)

    return

end


x = ...
    x(:).';


if isempty(t)

    error([ ...
        '%s is time-dependent, so a timestep index must be supplied ' ...
        'to transpirationDemand.'], ...
        name);

end


if ...
        ~isscalar(t) ...
        || t < 1 ...
        || t > numel(x) ...
        || t ~= round(t)

    error( ...
        'Invalid timestep index for %s.', ...
        name);

end


x = ...
    x(t);

end


%% =========================================================================
function es = saturationVaporPressureKPa(TC)

es = ...
    0.61078 ...
    .* exp( ...
        (17.2694 .* TC) ...
        ./ (TC + 237.29));

end