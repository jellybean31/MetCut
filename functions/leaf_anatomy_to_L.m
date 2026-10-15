function A = leaf_anatomy_to_L(anatomy)
%LEAF_ANATOMY_TO_L
%
% Derive the characteristic liquid-water transport geometry used in
% the leaf-water Peclet model.
%
% This function separates:
%
%   1. physical transport distance
%   2. tortuosity
%   3. active liquid-water flow area
%   4. fraction of transpiration routed through the pathway
%
%
% The Peclet number can then be written:
%
%   Pe = vPath * lPath / D
%
% where:
%
%   vPath =
%       E * flowFraction ...
%       / (Cwater * effectiveFlowAreaFraction)
%
%
% Equivalently:
%
%   Pe =
%       E * Leff / (Cwater * D)
%
% with:
%
%   Leff =
%       lPath ...
%       * flowFraction ...
%       / effectiveFlowAreaFraction
%
%
% REQUIRED
%
% anatomy.interveinalDistance_m
% anatomy.mesophyllThickness_m
% anatomy.effectiveFlowAreaFraction
% anatomy.flowFraction
% anatomy.tortuosity
%
%
% OPTIONAL
%
% anatomy.lateralPathFraction
%     Fraction of interveinal spacing defining the characteristic
%     lateral source-to-evaporation distance.
%
%     Default = 0.5
%
% anatomy.verticalPathFraction
%     Fraction of mesophyll thickness defining the characteristic
%     through-thickness transport distance.
%
%     Default = 0.5



%% ========================================================================
% Required inputs
% ========================================================================

required = { ...
    'interveinalDistance_m', ...
    'mesophyllThickness_m', ...
    'effectiveFlowAreaFraction', ...
    'flowFraction', ...
    'tortuosity'};


for k = 1:numel(required)

    if ~isfield(anatomy,required{k})

        error( ...
            'anatomy.%s is required.', ...
            required{k});

    end

end


dIV = ...
    anatomy.interveinalDistance_m;


tMes = ...
    anatomy.mesophyllThickness_m;


aEff = ...
    anatomy.effectiveFlowAreaFraction;


fFlow = ...
    anatomy.flowFraction;


tau = ...
    anatomy.tortuosity;


%% ========================================================================
% Optional geometric fractions
% ========================================================================

if isfield(anatomy,'lateralPathFraction')

    fLat = ...
        anatomy.lateralPathFraction;

else

    % Characteristic distance from vein to midpoint between veins.
    fLat = ...
        0.5;

end


if isfield(anatomy,'verticalPathFraction')

    fVert = ...
        anatomy.verticalPathFraction;

else

    % Development assumption:
    % distributed evaporation gives a characteristic distance
    % of approximately half the mesophyll thickness.
    fVert = ...
        0.5;

end


%% ========================================================================
% Validation
% ========================================================================

if ...
        ~isscalar(dIV) ...
        || ~isfinite(dIV) ...
        || dIV <= 0

    error( ...
        'interveinalDistance_m must be positive.');

end


if ...
        ~isscalar(tMes) ...
        || ~isfinite(tMes) ...
        || tMes <= 0

    error( ...
        'mesophyllThickness_m must be positive.');

end


if ...
        ~isscalar(aEff) ...
        || ~isfinite(aEff) ...
        || aEff <= 0 ...
        || aEff > 1

    error([ ...
        'effectiveFlowAreaFraction must lie in (0,1].']);

end


if ...
        ~isscalar(fFlow) ...
        || ~isfinite(fFlow) ...
        || fFlow <= 0 ...
        || fFlow > 1

    error( ...
        'flowFraction must lie in (0,1].');

end


if ...
        ~isscalar(tau) ...
        || ~isfinite(tau) ...
        || tau < 1

    error( ...
        'tortuosity must be >= 1.');

end


if ...
        ~isscalar(fLat) ...
        || ~isfinite(fLat) ...
        || fLat < 0 ...
        || fLat > 1

    error( ...
        'lateralPathFraction must lie in [0,1].');

end


if ...
        ~isscalar(fVert) ...
        || ~isfinite(fVert) ...
        || fVert < 0 ...
        || fVert > 1

    error( ...
        'verticalPathFraction must lie in [0,1].');

end


%% ========================================================================
% Physical source-to-evaporation geometry
% ========================================================================

lateralDistance = ...
    fLat ...
    .* dIV;


verticalDistance = ...
    fVert ...
    .* tMes;


directPathLength = ...
    hypot( ...
        lateralDistance, ...
        verticalDistance);


%% ========================================================================
% Tortuous liquid-water path
% ========================================================================

pathLength = ...
    tau ...
    .* directPathLength;


%% ========================================================================
% Flow-concentration scaling
%
% E is expressed per unit projected leaf area.
%
% If liquid flow is concentrated into only a fraction aEff of that
% area, the pathway-scale velocity increases by 1/aEff.
%
% flowFraction allows only part of total transpiration to be assigned
% to this pathway.
% ========================================================================

velocityScale = ...
    fFlow ...
    ./ aEff;


%% ========================================================================
% Backward-compatible effective length
% ========================================================================

Leff = ...
    pathLength ...
    .* velocityScale;


%% ========================================================================
% Outputs
% ========================================================================

A = ...
    struct;


A.interveinalDistance_m = ...
    dIV;


A.mesophyllThickness_m = ...
    tMes;


A.lateralPathFraction = ...
    fLat;


A.verticalPathFraction = ...
    fVert;


A.lateralDistance_m = ...
    lateralDistance;


A.verticalDistance_m = ...
    verticalDistance;


A.directPathLength_m = ...
    directPathLength;


A.tortuosity = ...
    tau;


A.pathLength_m = ...
    pathLength;


A.effectiveFlowAreaFraction = ...
    aEff;


A.flowFraction = ...
    fFlow;


A.velocityScale = ...
    velocityScale;


% Retained so existing code can still use:
%
%   Pe = E * Leff / (C * D)

A.Leff_m = ...
    Leff;

end