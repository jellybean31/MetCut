function B = effective7_structural_basis()
%EFFECTIVE7_STRUCTURAL_BASIS
% Load the architecture-v2 raw rank-4 structural basis.
%
% The v2 basis spans the four-dimensional raw observable subspace in:
%   E01 E02 E03 E04 E05 stage6 stage7
%
% All seven coordinates have per-mil units.
%
% Expected package reference:
%   data/reference/effective7_architecture_v2_structural_basis_v1.mat


thisFile=mfilename('fullpath');
functionDir=fileparts(thisFile);
projectRoot=fileparts(functionDir);
checkpointFile=fullfile( ...
    projectRoot,'data','reference', ...
    'effective7_architecture_v2_structural_basis_v1.mat');

if ~isfile(checkpointFile)
    error(['Accepted architecture-v2 structural-basis reference not found:\n%s'], ...
           checkpointFile);
end

S=load(checkpointFile,'C');

if ~isfield(S,'C')
    error('Checkpoint does not contain expected structure C.');
end

C=S.C;

expectedVersion='2026-08-29_effective7_architecture_v2_structural_basis_v1';
if ~isfield(C,'version') || ~strcmp(C.version,expectedVersion)
    error('Unexpected architecture-v2 structural-basis checkpoint version.');
end

if ~isfield(C,'validation') || ...
        ~isfield(C.validation,'acceptance') || ...
        ~C.validation.acceptance.overall
    error('architecture-v2 basis checkpoint is not accepted.');
end

B=struct;
B.version=C.version;
B.architectureVersion=C.architectureVersion;
B.parameterIDs=C.parameterIDs;
B.parameterNames=C.parameterNames;
B.coordinateUnits=C.coordinateUnits;
B.referenceEffective=C.referenceEffective7;
B.dimension=C.structuralDimension;
B.V4=C.rawStructuralBasis;
B.projector=C.rawStructuralProjector;
B.formula=C.formula;
B.modeLoadingTable=C.modeLoadingTable;
B.priorStatus=C.priorStatus;
B.checkpointFile=checkpointFile;

end
