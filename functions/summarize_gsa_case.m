function O = summarize_gsa_case(M,B)
%SUMMARIZE_GSA_CASE Compact scalar outputs for sensitivity analysis.
%
% Diagnostics:
%   max/rms |cytosol - intracellular source|
%   max/rms |plastid - intracellular source|
%   max/rms |plastid - cytosol|


O=struct;

O.finalRootSource_permil=M.rootBoundary.deltaDNode_permil(end);
O.finalXylemDelivered_permil=M.xylemIsotopes.leafBoundary.deltaDWater(end);
O.finalLeafSymplast_permil=M.leafIsotopes.deltaDSymplastBulk(end);
O.finalEvaporationSite_permil=M.leafIsotopes.deltaDEvaporationSite(end);
O.finalCytosol_permil=M.intracellular.deltaD.cytosol(end);
O.finalPlastid_permil=M.intracellular.deltaD.plastid(end);

A=M.biochemistry.alkaneCarbonNumber(:);
W=M.biochemistry.wax.deltaDMatrix;

for C=[27 29 31 33]
    ii=find(A==C,1);
    if isempty(ii), continue; end
    O.(sprintf('finalC%d_permil',C))=W(ii,end);
end

if isfield(O,'finalC33_permil') && isfield(O,'finalC27_permil')
    O.finalC33minusC27_permil=O.finalC33_permil-O.finalC27_permil;
end

if isfield(O,'finalC31_permil') && isfield(O,'finalC29_permil')
    O.finalC31minusC29_permil=O.finalC31_permil-O.finalC29_permil;
end

% -------------------------------------------------------------------------
% Transient intracellular memory diagnostics
% -------------------------------------------------------------------------

cy=M.intracellular.deltaD.cytosol(:);
pl=M.intracellular.deltaD.plastid(:);

if isfield(M.intracellular,'sourceDeltaDWater')
    src=M.intracellular.sourceDeltaDWater(:);
else
    src=M.leafIsotopes.deltaDSymplastBulk(:);
end

if numel(src)~=numel(cy)
    error('Intracellular source and cytosol trajectories have different lengths.');
end

dCy=cy-src;
dPl=pl-src;
dPlCy=pl-cy;

O.maxAbsCytosolMinusSource_permil=max(abs(dCy),[],'omitnan');
O.rmsCytosolMinusSource_permil=sqrt(mean(dCy.^2,'omitnan'));

O.maxAbsPlastidMinusSource_permil=max(abs(dPl),[],'omitnan');
O.rmsPlastidMinusSource_permil=sqrt(mean(dPl.^2,'omitnan'));

O.maxAbsPlastidMinusCytosol_permil=max(abs(dPlCy),[],'omitnan');
O.rmsPlastidMinusCytosol_permil=sqrt(mean(dPlCy.^2,'omitnan'));

% -------------------------------------------------------------------------
% Existing water-age diagnostics
% -------------------------------------------------------------------------

if isfield(M.xylemIsotopes,'leafBoundary') && ...
        isfield(M.xylemIsotopes.leafBoundary,'meanWaterAgeDays')

    age=M.xylemIsotopes.leafBoundary.meanWaterAgeDays(:);

    O.finalXylemWaterAge_days=age(end);
    O.maxXylemWaterAge_days=max(age,[],'omitnan');
    O.meanXylemWaterAge_days=mean(age,'omitnan');
end

if isfield(M.leafIsotopes,'pecletTransferFactorSpatialObserved')
    phi=M.leafIsotopes.pecletTransferFactorSpatialObserved;
    O.maxLeafPhi=max(phi,[],'all','omitnan');
    O.finalLeafPhi=phi(end);
end

% -------------------------------------------------------------------------
% Closure diagnostics
% -------------------------------------------------------------------------

O.maxXylemDClosure_molD_m2=max(abs( ...
    M.xylemIsotopes.deuteriumInventoryClosureResidual_molD_m2), ...
    [],'all','omitnan');

O.maxLeafDClosure_molD_m2=max(abs( ...
    M.leafIsotopes.DInventoryClosureResidual_molD_m2), ...
    [],'all','omitnan');

O.maxIntracellularDClosure_molD_m2=max(abs( ...
    M.intracellular.wholeSystemDClosureResidual_molD_m2), ...
    [],'all','omitnan');

end
