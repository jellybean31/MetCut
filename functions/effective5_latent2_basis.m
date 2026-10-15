
function B = effective5_latent2_basis(Truth)
%EFFECTIVE5_LATENT2_BASIS
% Define reusable two-dimensional biochemical latent coordinates.
%
% Reference:
%   midpoint of the admissible E01-E05 effective-state ranges.
%
% Scaling:
%   z_j = (E_j - Eref_j) / range_j
%
% Coordinates:
%   L1 = v1' * z
%   L2 = v2' * z
%
% where v1 and v2 are the first two right-singular vectors of the scaled
% 4x5 homolog-sensitivity Jacobian evaluated at Eref.


D=effective5_parameter_bounds(Truth.model.intracellular);

lo=[D.low].';
hi=[D.high].';
range=hi-lo;
eref=0.5.*(lo+hi);

R=effective5_practical_rank_subspace_at_point( ...
    Truth,eref,'latent2_reference');

V=R.rightSingularVectors;
V2=V(:,1:2);

% Deterministic signs matching the interpretation:
% mode 1 -> positive common shift
% mode 2 -> positive ER contribution if possible.
%
% Mode 1:
if sum(V2(:,1))<0
    V2(:,1)=-V2(:,1);
end
% Mode 2:
if V2(5,2)<0
    V2(:,2)=-V2(:,2);
end

B=struct;
B.version='2026-08-28_effective5_latent2_basis_v1';
B.parameterDefinitions=D;
B.parameterIDs=string({D.id}).';
B.parameterNames=string({D.name}).';
B.referenceEffective=eref;
B.parameterRange=range;
B.V2=V2;
B.fullRightSingularVectors=V;
B.referenceAnalysis=R;

B.formula={ ...
    'z = (E - Eref) ./ range'; ...
    'L1 = v1'' * z'; ...
    'L2 = v2'' * z'; ...
    'E_representative = Eref + range .* (v1*L1 + v2*L2)'};

B.mode1Table=table( ...
    B.parameterIDs,B.parameterNames,V2(:,1), ...
    'VariableNames',{'ParameterID','Parameter','L1Loading_Standardized'});

B.mode2Table=table( ...
    B.parameterIDs,B.parameterNames,V2(:,2), ...
    'VariableNames',{'ParameterID','Parameter','L2Loading_Standardized'});

fprintf('\n');
fprintf('============================================================\n');
fprintf('EFFECTIVE-5 LATENT-2 BASIS\n');
fprintf('============================================================\n');
fprintf('Reference = midpoint of admissible ranges\n');

fprintf('\nL1 standardized loadings\n');
disp(B.mode1Table);

fprintf('\nL2 standardized loadings\n');
disp(B.mode2Table);

fprintf('============================================================\n');

end
