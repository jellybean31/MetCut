function R = reconstruct_richards_interval_transfers(S)
%RECONSTRUCT_RICHARDS_INTERVAL_TRANSFERS
%
% Reconstruct the NET water transfer across every internal Richards
% cell face over each stored coarse interval using the finite-volume water
% balance already contained in S.
%
% Sign conventions:
%
%   interfaceTransferNet_mol_m2 > 0 : downward transfer i -> i+1
%   interfaceTransferNet_mol_m2 < 0 : upward transfer   i+1 -> i
%
%   rootTransferSigned_mol_m2 > 0   : soil -> root uptake
%   rootTransferSigned_mol_m2 < 0   : root -> soil release


%% ========================================================================
% Required fields
% ========================================================================

required = { ...
    'timeDays', ...
    'waterAmount_mol_m2', ...
    'JinfiltrationActualInterval_mol_m2_s', ...
    'Jdeep_mol_m2_s', ...
    'JrootSignedInterval_mol_m2_s'};

for q = 1:numel(required)

    if ~isfield(S,required{q})

        error( ...
            'Missing Richards output field: %s', ...
            required{q});

    end

end


%% ========================================================================
% Dimensions
% ========================================================================

t = ...
    S.timeDays(:).';

nT = ...
    numel(t);

nInt = ...
    nT - 1;

W = ...
    S.waterAmount_mol_m2;

[nCell,nTW] = ...
    size(W);

if nTW ~= nT

    error( ...
        'waterAmount_mol_m2 dimensions are inconsistent with timeDays.');

end


dt_s = ...
    diff(t) ...
    .* 86400;


%% ========================================================================
% Integrated coarse-interval terms
% ========================================================================

infiltration = ...
    reshape( ...
        S.JinfiltrationActualInterval_mol_m2_s, ...
        1,[]) ...
    .* dt_s;


deepRecorded = ...
    reshape( ...
        S.Jdeep_mol_m2_s, ...
        1,[]) ...
    .* dt_s;


Jroot = ...
    S.JrootSignedInterval_mol_m2_s;

if ~isequal(size(Jroot),[nCell nInt])

    error( ...
        'JrootSignedInterval_mol_m2_s dimensions are inconsistent.');

end


rootSigned = ...
    Jroot ...
    .* reshape(dt_s,1,[]);


dW = ...
    diff(W,1,2);


%% ========================================================================
% Internal-face reconstruction
% ========================================================================

if nCell > 1

    faceNet = ...
        zeros(nCell-1,nInt);


    % Top cell:
    %
    %   dW1 = infiltration - F1 - root1

    faceNet(1,:) = ...
        infiltration ...
        - rootSigned(1,:) ...
        - dW(1,:);


    % Interior recurrence:
    %
    %   dWi = F(i-1) - Fi - root_i

    for i = 2:nCell-1

        faceNet(i,:) = ...
            faceNet(i-1,:) ...
            - rootSigned(i,:) ...
            - dW(i,:);

    end


    % Bottom cell:
    %
    %   dWN = F(N-1) - deep - root_N

    deepReconstructed = ...
        faceNet(end,:) ...
        - rootSigned(end,:) ...
        - dW(end,:);

else

    faceNet = ...
        zeros(0,nInt);


    deepReconstructed = ...
        infiltration ...
        - rootSigned(1,:) ...
        - dW(1,:);

end


%% ========================================================================
% Cellwise reconstruction residual
% ========================================================================

cellResidual = ...
    zeros(nCell,nInt);


if nCell == 1

    predicted = ...
        infiltration ...
        - deepReconstructed ...
        - rootSigned(1,:);


    cellResidual(1,:) = ...
        dW(1,:) ...
        - predicted;

else

    predictedTop = ...
        infiltration ...
        - faceNet(1,:) ...
        - rootSigned(1,:);


    cellResidual(1,:) = ...
        dW(1,:) ...
        - predictedTop;


    for i = 2:nCell-1

        predictedInterior = ...
            faceNet(i-1,:) ...
            - faceNet(i,:) ...
            - rootSigned(i,:);


        cellResidual(i,:) = ...
            dW(i,:) ...
            - predictedInterior;

    end


    predictedBottom = ...
        faceNet(end,:) ...
        - deepReconstructed ...
        - rootSigned(end,:);


    cellResidual(end,:) = ...
        dW(end,:) ...
        - predictedBottom;

end


%% ========================================================================
% Compare reconstructed lower boundary with stored Richards diagnostic
% ========================================================================

deepDifference = ...
    deepReconstructed ...
    - deepRecorded;


if isfield(S,'stepMassBalanceResidual_mol_m2')

    expectedDeepDifference = ...
        -reshape( ...
            S.stepMassBalanceResidual_mol_m2, ...
            1,[]);


    deepDifferenceConsistency = ...
        deepDifference ...
        - expectedDeepDifference;

else

    expectedDeepDifference = ...
        NaN(1,nInt);


    deepDifferenceConsistency = ...
        NaN(1,nInt);

end


%% ========================================================================
% Outputs
% ========================================================================

R.timeDays = ...
    t;

R.dt_s = ...
    dt_s;

R.infiltration_mol_m2 = ...
    infiltration;

R.rootTransferSigned_mol_m2 = ...
    rootSigned;

R.interfaceTransferNet_mol_m2 = ...
    faceNet;

R.deepRecorded_mol_m2 = ...
    deepRecorded;

R.deepReconstructed_mol_m2 = ...
    deepReconstructed;

R.deepDifference_mol_m2 = ...
    deepDifference;

R.expectedDeepDifferenceFromMassResidual_mol_m2 = ...
    expectedDeepDifference;

R.deepDifferenceConsistency_mol_m2 = ...
    deepDifferenceConsistency;

R.cellReconstructionResidual_mol_m2 = ...
    cellResidual;

R.maxAbsCellReconstructionResidual_mol_m2 = ...
    max(abs(cellResidual),[],'all');

R.maxAbsDeepDifference_mol_m2 = ...
    max(abs(deepDifference),[],'all');

R.maxAbsDeepDifferenceConsistency_mol_m2 = ...
    max(abs(deepDifferenceConsistency),[],'all','omitnan');

end
