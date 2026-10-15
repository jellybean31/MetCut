function Q = construct_rooting_arch_hydraulic_reference_summary(Mreference)
%CONSTRUCT_ROOTING_ARCH_HYDRAULIC_REFERENCE_SUMMARY_CANDIDATE
%
% Candidate Helper B implementation.
%
%   Q = construct_rooting_arch_hydraulic_reference_summary(Mreference)
%
% Summary constructor for an already-supplied wrapped hydraulic reference.
%
% Sources:
%   Mreference.hydraulicCheckpointPath
%   Mreference.grid.rootLengthGridByZone_m_m2
%   Mreference.grid.originalZoneBottom_m
%
% Diagnostic structural cross-check:
%   Mreference.grid.rootLengthOriginal_m_m2
%
% Structural fractions:
%   rootLengthGridByZone_m_m2 / sum(rootLengthGridByZone_m_m2)
%
% Zone bounds:
%   [0 z1; z1 z2; z2 z3]

Q=struct;
Q.version="2026-09-06_rooting_helperB_candidate_v1";
Q.status="indeterminate";
Q.hydraulicCheckpointPath="";
Q.zoneUptakeFractions=[NaN NaN NaN];
Q.zoneDepthBounds_m=nan(3,2);

Q.diagnostics=struct;
Q.diagnostics.reason="";
Q.diagnostics.structuralSource="Mreference.grid.rootLengthGridByZone_m_m2";
Q.diagnostics.structuralDiagnosticSource="Mreference.grid.rootLengthOriginal_m_m2";
Q.diagnostics.geometrySource="Mreference.grid.originalZoneBottom_m";
Q.diagnostics.fractionClosureResidual=NaN;
Q.diagnostics.gridVsOriginalFractionResidual=NaN;
Q.diagnostics.structuralFractionFinite=false;
Q.diagnostics.geometryFinite=false;

if ~isstruct(Mreference) || ~isscalar(Mreference)
    Q.diagnostics.reason="invalid_Mreference";
    return
end

%% Hydraulic reference identity
if ~isfield(Mreference,'hydraulicCheckpointPath')
    Q.diagnostics.reason="missing_hydraulicCheckpointPath";
    return
end

id=string(Mreference.hydraulicCheckpointPath);

if ~isscalar(id) || strlength(id)==0
    Q.diagnostics.reason="missing_hydraulicCheckpointPath";
    return
end

Q.hydraulicCheckpointPath=id;

%% Structural rooting source
if ~isfield(Mreference,'grid') || ...
        ~isstruct(Mreference.grid) || ...
        ~isscalar(Mreference.grid)
    Q.diagnostics.reason="missing_grid";
    return
end

if ~isfield(Mreference.grid,'rootLengthGridByZone_m_m2')
    Q.diagnostics.reason="missing_rootLengthGridByZone_m_m2";
    return
end

L=Mreference.grid.rootLengthGridByZone_m_m2;

if ~isnumeric(L) || numel(L)~=3
    Q.diagnostics.reason="invalid_rootLengthGridByZone_shape";
    return
end

L=double(reshape(L,1,3));

if any(~isfinite(L)) || any(L<0) || sum(L)<=0
    Q.diagnostics.reason="invalid_rootLengthGridByZone_values";
    return
end

f=L./sum(L);
closure=abs(sum(f)-1);

Q.diagnostics.fractionClosureResidual=closure;
Q.diagnostics.structuralFractionFinite=all(isfinite(f));

if closure>64*eps(max([1 abs(f)]))
    Q.diagnostics.reason="structural_fraction_nonclosure";
    return
end

%% Optional diagnostic cross-check against original zone root lengths
if isfield(Mreference.grid,'rootLengthOriginal_m_m2')
    L0=Mreference.grid.rootLengthOriginal_m_m2;

    if isnumeric(L0) && numel(L0)==3
        L0=double(reshape(L0,1,3));

        if all(isfinite(L0)) && all(L0>=0) && sum(L0)>0
            f0=L0./sum(L0);
            Q.diagnostics.gridVsOriginalFractionResidual=max(abs(f-f0));
        end
    end
end

%% Zone geometry
if ~isfield(Mreference.grid,'originalZoneBottom_m')
    Q.diagnostics.reason="missing_originalZoneBottom_m";
    return
end

zb=Mreference.grid.originalZoneBottom_m;

if ~isnumeric(zb) || numel(zb)~=3
    Q.diagnostics.reason="invalid_originalZoneBottom_shape";
    return
end

zb=double(reshape(zb,1,3));

if any(~isfinite(zb)) || any(zb<=0) || any(diff(zb)<=0)
    Q.diagnostics.reason="invalid_originalZoneBottom_values";
    return
end

bounds=[
    0     zb(1)
    zb(1) zb(2)
    zb(2) zb(3)
];

Q.diagnostics.geometryFinite=all(isfinite(bounds(:)));

%% Final resolved summary
Q.status="resolved";
Q.zoneUptakeFractions=f;
Q.zoneDepthBounds_m=bounds;
Q.diagnostics.reason="ok";

end
