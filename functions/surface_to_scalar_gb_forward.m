function out = surface_to_scalar_gb_forward(gsScalar,qSurface,gbSurface)
% SURFACE_TO_SCALAR_GB_FORWARD
%
% Forward aggregation from scalar hydraulic stomatal conductance plus
% surface allocation and surface boundary conductance to a flux-equivalent
% scalar boundary conductance.
%
% Surface order:
%   [ABAXIAL, ADAXIAL]
%
% INPUTS
% ------
% gsScalar
%   finite nonnegative scalar, mol m^-2 s^-1
%   reference area: ONE_SIDED_LAMINA_PLANFORM_AREA
%
% qSurface
%   finite nonnegative two-element scale-free allocation vector
%   physiology/provenance must be authorized by the caller
%
% gbSurface
%   finite nonnegative two-element boundary conductance vector,
%   mol m^-2 s^-1
%   reference area: ONE_SIDED_LAMINA_PLANFORM_AREA
%
% OUTPUT
% ------
% out.surfaceOrder
% out.qHat
% out.gsScalar
% out.gsSurface
% out.gbSurface
% out.gwSurface
% out.gwTotal
% out.gbEquivalent
% out.equivalentExists
% out.equivalentIdentifiable
% out.referenceArea
% out.units
%
% NUMERICAL NOTE
% --------------
% The direct inverse:
%
%   gbEquivalent = gsScalar*gwTotal/(gsScalar-gwTotal)
%
% is algebraically correct but suffers subtractive cancellation when
% gwTotal/gsScalar -> 1.


assert(isnumeric(gsScalar) && isscalar(gsScalar) && ...
    isfinite(gsScalar) && gsScalar>=0, ...
    'gsScalar must be a finite nonnegative numeric scalar.');

assert(isnumeric(qSurface) && numel(qSurface)==2, ...
    'qSurface must be a numeric two-element vector.');

qSurface=double(qSurface(:).');

assert(all(isfinite(qSurface)) && all(qSurface>=0) && sum(qSurface)>0, ...
    'qSurface must be finite, nonnegative, and have positive sum.');

assert(isnumeric(gbSurface) && numel(gbSurface)==2, ...
    'gbSurface must be a numeric two-element vector.');

gbSurface=double(gbSurface(:).');

assert(all(isfinite(gbSurface)) && all(gbSurface>=0), ...
    'gbSurface must be finite and nonnegative.');

qHat=qSurface./sum(qSurface);
gsSurface=gsScalar.*qHat;

gwSurface=zeros(1,2);

for s=1:2
    if gsSurface(s)==0 || gbSurface(s)==0
        gwSurface(s)=0;
    else
        gwSurface(s)= ...
            gsSurface(s).*gbSurface(s)./ ...
            (gsSurface(s)+gbSurface(s));
    end
end

gwTotal=sum(gwSurface);

out=struct();
out.surfaceOrder=["ABAXIAL","ADAXIAL"];
out.qHat=qHat;
out.gsScalar=gsScalar;
out.gsSurface=gsSurface;
out.gbSurface=gbSurface;
out.gwSurface=gwSurface;
out.gwTotal=gwTotal;
out.referenceArea="ONE_SIDED_LAMINA_PLANFORM_AREA";
out.units="mol m^-2 s^-1";

if gsScalar==0
    assert(gwTotal==0, ...
        'Internal error: zero gsScalar produced nonzero gwTotal.');

    out.gbEquivalent=NaN;
    out.equivalentExists=true;
    out.equivalentIdentifiable=false;
    return
end

activeSurface = qHat>0;
if nnz(activeSurface)==1
    iActive=find(activeSurface,1);
    out.gbEquivalent=gbSurface(iActive);
    out.equivalentExists=true;
    out.equivalentIdentifiable=true;
    return
end

gsSurfaceStable = gsScalar.*qHat;
stableEquivalentGap = sum( ...
    gsSurfaceStable.^2 ./ (gsSurfaceStable + gbSurface));
out.equivalentExists = ...
    isfinite(gwTotal) && gwTotal>=0 && ...
    isfinite(stableEquivalentGap) && stableEquivalentGap>0;

out.equivalentIdentifiable = ...
    out.equivalentExists;

if ~out.equivalentExists
    out.gbEquivalent=NaN;
    return
end

if gwTotal==0
    out.gbEquivalent=0;
    return
end

active=find(qHat>0);

if numel(active)==1
    out.gbEquivalent=gbSurface(active);
    return
end

den=gsScalar.*qHat + gbSurface;
term=zeros(1,2);

for s=1:2
    if qHat(s)==0
        term(s)=0;
    else
        assert(den(s)>0, ...
            'Internal error: active surface has nonpositive denominator.');
        term(s)=qHat(s).^2./den(s);
    end
end

S=sum(term);
r=gwTotal./gsScalar;

assert(isfinite(S) && S>0, ...
    'Internal error: equivalent-conductance denominator is invalid.');

out.gbEquivalent=r./S;

end
