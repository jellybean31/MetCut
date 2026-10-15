
function R = effective5_practical_rank_subspace_at_point(Truth,x,label,options)
%EFFECTIVE5_PRACTICAL_RANK_SUBSPACE_AT_POINT
% Practical-rank and identifiable-subspace analysis for the five-state
% biochemical inverse.
%
% Builds the same dimensionless scaled Jacobian used previously:
%
%   G = diag(1/sigmaWax) * J * diag(parameterRange)
%
% and analyzes:
%   - all singular values at full precision
%   - singular-value ratios
%   - exact numerical rank
%   - practical rank under user-defined relative thresholds
%   - rank-k reconstruction error for k=1..4
%   - dominant identifiable right-singular modes
%   - near-null subspace basis
%
% options fields:
%   sigmaWax_permil          default Truth.config.biochemistry.sigmaWax_permil
%   relativeRankThresholds   default [1e-1 1e-2 1e-3 1e-6 1e-9]
%   finiteDifferenceStep     default 1e-4 of parameter range

if nargin<3 || isempty(label), label='point'; end
if nargin<4 || isempty(options), options=struct; end

if ~isfield(options,'sigmaWax_permil')
    options.sigmaWax_permil=Truth.config.biochemistry.sigmaWax_permil;
end
if ~isfield(options,'relativeRankThresholds')
    options.relativeRankThresholds=[1e-1 1e-2 1e-3 1e-6 1e-9];
end
if ~isfield(options,'finiteDifferenceStep')
    options.finiteDifferenceStep=1e-4;
end

fd=struct;
fd.stepFractionOfRange=options.finiteDifferenceStep;
fd.minimumStep_permil=1e-5;

J=effective5_finite_difference_jacobian(Truth,x,fd);

D=effective5_parameter_bounds(Truth.model.intracellular);
ranges=([D.high]-[D.low]).';

sigma=options.sigmaWax_permil;
if isscalar(sigma)
    sigma=repmat(sigma,size(J,1),1);
else
    sigma=sigma(:);
end

G=diag(1./sigma)*J*diag(ranges);

[U,S,V]=svd(G,'econ');
s=diag(S);

% For full 4x5 right-singular basis, obtain null complement too.
[Uf,Sf,Vf]=svd(G); 
sFull=diag(Sf(:,1:size(G,1)));

sRel=s./s(1);
adjacentRatio=NaN(size(s));
adjacentRatio(1)=1;
for k=2:numel(s)
    adjacentRatio(k)=s(k)/s(k-1);
end

% Numerical rank.
tol=max(size(G))*eps(max(s))*max(s);
numericalRank=sum(s>tol);

% Practical rank under several explicit relative thresholds.
thr=options.relativeRankThresholds(:);
prank=NaN(numel(thr),1);
for k=1:numel(thr)
    prank(k)=sum(sRel>=thr(k));
end
PracticalRankTable=table(thr,prank, ...
    'VariableNames',{'RelativeThreshold','PracticalRank'});

% Rank-k reconstruction diagnostics.
kmax=numel(s);
RankK=(1:kmax).';
RelativeFrobeniusError=NaN(kmax,1);
ExplainedSquaredSensitivity=NaN(kmax,1);
MaxAbsScaledJacobianError=NaN(kmax,1);

ssTotal=sum(s.^2);

for k=1:kmax
    Gk=U(:,1:k)*S(1:k,1:k)*V(:,1:k).';
    d=G-Gk;
    RelativeFrobeniusError(k)=norm(d,'fro')/norm(G,'fro');
    ExplainedSquaredSensitivity(k)=sum(s(1:k).^2)/ssTotal;
    MaxAbsScaledJacobianError(k)=max(abs(d),[],'all');
end

ReconstructionTable=table( ...
    RankK,RelativeFrobeniusError,ExplainedSquaredSensitivity, ...
    MaxAbsScaledJacobianError);

% Mode tables: standardized coordinates use V directly.
paramIDs=string({D.id}).';
paramNames=string({D.name}).';

modeCount=size(V,2);
ModeTables=cell(modeCount,1);

for k=1:modeCount
    vStd=V(:,k);

    % Convert dz to physical dx = range .* dz, then normalize for display.
    vPhys=ranges.*vStd;
    vPhys=vPhys./norm(vPhys);

    % Deterministic sign for readability.
    [~,imax]=max(abs(vStd));
    if vStd(imax)<0
        vStd=-vStd;
        vPhys=-vPhys;
    end

    ModeTables{k}=table( ...
        paramIDs,paramNames,vStd,vPhys, ...
        'VariableNames',{ ...
        'ParameterID','Parameter', ...
        'Loading_Standardized', ...
        'Loading_PhysicalNormalized'});
end

% Full right-singular basis contains fifth exact-null vector.
vNull=Vf(:,end);
vNullPhys=ranges.*vNull;
vNullPhys=vNullPhys./norm(vNullPhys);
[~,imax]=max(abs(vNull));
if vNull(imax)<0
    vNull=-vNull;
    vNullPhys=-vNullPhys;
end

ExactNullTable=table( ...
    paramIDs,paramNames,vNull,vNullPhys, ...
    'VariableNames',{ ...
    'ParameterID','Parameter', ...
    'Loading_Standardized', ...
    'Loading_PhysicalNormalized'});

% Near-null subspace candidates: report the bottom three right-singular
% directions of the full V basis. If practical rank is ~2, these span the
% approximately unidentifiable subspace.
NearNullBasisStd=Vf(:,3:5);
NearNullBasisPhys=diag(ranges)*NearNullBasisStd;

for k=1:size(NearNullBasisPhys,2)
    NearNullBasisPhys(:,k)= ...
        NearNullBasisPhys(:,k)./norm(NearNullBasisPhys(:,k));
end

R=struct;
R.version='2026-08-28_effective5_practical_rank_subspace_v1';
R.label=label;
R.x=x(:);
R.rawJacobian=J;
R.scaledJacobian=G;
R.singularValues=s;
R.relativeSingularValues=sRel;
R.adjacentSingularRatios=adjacentRatio;
R.numericalRank=numericalRank;
R.numericalTolerance=tol;
R.practicalRankTable=PracticalRankTable;
R.reconstructionTable=ReconstructionTable;
R.rightSingularVectors=Vf;
R.modeTables=ModeTables;
R.exactNullTable=ExactNullTable;
R.nearNullBasisStandardized=NearNullBasisStd;
R.nearNullBasisPhysical=NearNullBasisPhys;
R.parameterRanges=ranges;
R.sigmaWax_permil=sigma;

fprintf('\n');
fprintf('============================================================\n');
fprintf('EFFECTIVE-5 PRACTICAL RANK / SUBSPACE: %s\n',upper(label));
fprintf('============================================================\n');

fprintf('Singular values (full precision):\n');
for k=1:numel(s)
    fprintf('  s%d = %.16e | relative %.16e | adjacent ratio %.16e\n', ...
        k,s(k),sRel(k),adjacentRatio(k));
end

fprintf('Numerical rank: %d / 4\n',numericalRank);

fprintf('\nPractical rank under explicit relative thresholds\n');
disp(PracticalRankTable);

fprintf('\nRank-k reconstruction\n');
disp(ReconstructionTable);

for k=1:min(2,numel(ModeTables))
    fprintf('\nIdentifiable mode %d | singular value %.16e\n',k,s(k));
    disp(ModeTables{k});
end

fprintf('\nExact algebraic null direction\n');
disp(ExactNullTable);

fprintf('============================================================\n');

end
