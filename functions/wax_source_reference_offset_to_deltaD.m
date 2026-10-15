function [deltaSource_permil,alphaSource] = ...
    wax_source_reference_offset_to_deltaD( ...
    referenceOffset_permil,localWater_permil,referenceWater_permil)
%WAX_SOURCE_REFERENCE_OFFSET_TO_DELTAD
% Convert a source-state prior expressed as a deltaD offset at a reference
% water isotope composition into an exact constant fractionation factor,
% then propagate it to arbitrary local-water deltaD.
%
% Existing biochemical source priors retain their numerical bounds as
% effective offsets at the standardized reference state:
%
%   Delta_ref = deltaD_source(ref) - deltaD_water(ref)
%
% They are converted internally to:
%
%   alpha_source =
%       (1 + (deltaD_ref + Delta_ref)/1000) /
%       (1 + deltaD_ref/1000)
%
% and propagated as:
%
%   1 + deltaD_source/1000 =
%       alpha_source * (1 + deltaD_localWater/1000)
%
% Therefore localWater == referenceWater exactly reproduces:
%
%   deltaD_source = referenceWater + referenceOffset
%
% Inputs may be scalars or arrays compatible by implicit expansion.

if nargin<3 || isempty(referenceWater_permil)
    referenceWater_permil=-100;
end

if ~isscalar(referenceWater_permil) || ~isfinite(referenceWater_permil) || ...
        referenceWater_permil<=-1000
    error('referenceWater_permil must be a finite scalar greater than -1000.');
end
if any(~isfinite(referenceOffset_permil(:))) || ...
        any(~isfinite(localWater_permil(:))) || ...
        any(localWater_permil(:)<=-1000)
    error('Offsets/water values must be finite; water deltaD must exceed -1000 per mil.');
end

sourceAtRef=referenceWater_permil + referenceOffset_permil;
if any(sourceAtRef(:)<=-1000)
    error('Reference offset produces a nonphysical source deltaD <= -1000 per mil.');
end

rRef=1+referenceWater_permil/1000;
alphaSource=(1+sourceAtRef/1000)./rRef;

if any(~isfinite(alphaSource(:))) || any(alphaSource(:)<=0)
    error('Computed source alpha must be positive and finite.');
end

deltaSource_permil=1000.*(alphaSource.*(1+localWater_permil/1000)-1);

if any(~isfinite(deltaSource_permil(:))) || any(deltaSource_permil(:)<=-1000)
    error('Propagated source deltaD must remain finite and greater than -1000 per mil.');
end
end
