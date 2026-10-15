function S = alkane_biosynthesis_h_stoichiometry_reaction_explicit(alkaneCarbonNumber,routeFractions)
%ALKANE_BIOSYNTHESIS_H_STOICHIOMETRY_REACTION_EXPLICIT
% Reaction-explicit H-atom bookkeeping for odd-chain plant n-alkanes.
%
% This is the reaction-resolution successor to
% alkane_biosynthesis_h_stoichiometry_postmalonate.m. It preserves the
% acetyl / malonyl / CER3 structural revision while separating
% reductive and protonation H by the reaction that installs them.
%
% Explicit H classes:
%   acetylStarterH
%   plastidRetainedMalonylH
%   erRetainedMalonylH
%   plastidKARNADPHHydrideH
%   plastidENRNADHHydrideH
%   plastidENRProtonH
%   erKCRNADPHHydrideH
%   erECRNADPHHydrideH
%   erECRProtonH
%   cer3AldehydeHydrideH
%
% Per plastid elongation cycle:
%   KAR hydride       = 1 H from plastid NADPH
%   ENR hydride       = 1 H from plastid NADH
%   ENR proton        = 1 H from plastid-local exchangeable/water H
%
% Per ER elongation cycle:
%   KCR hydride       = 1 H from ER/cytosolic NADPH
%   ECR hydride       = 1 H from ER/cytosolic NADPH
%   ECR proton        = 1 H from cytosol/ER-local exchangeable/water H
%
% The CER3 aldehydic H is retained into the alkane.

n = alkaneCarbonNumber(:).';
if isempty(n) || any(~isfinite(n)) || any(n~=round(n)) || ...
        any(mod(n,2)~=1) || any(n<17)
    error('alkaneCarbonNumber values must be odd integers >= 17.');
end

if nargin < 2 || isempty(routeFractions)
    f16 = 0; f18 = 1;
elseif isstruct(routeFractions)
    if ~isfield(routeFractions,'C16') || ~isfield(routeFractions,'C18')
        error('routeFractions structure must contain C16 and C18.');
    end
    f16 = routeFractions.C16; f18 = routeFractions.C18;
elseif isnumeric(routeFractions) && numel(routeFractions)==2
    f16 = routeFractions(1); f18 = routeFractions(2);
else
    error('routeFractions must be struct(C16,C18) or numeric [f16 f18].');
end

if ~isscalar(f16) || ~isscalar(f18) || ~isfinite(f16) || ~isfinite(f18) || ...
        f16<0 || f18<0 || abs(f16+f18-1)>1e-12
    error('Route fractions must be finite, nonnegative, and sum to 1.');
end

sourceClasses = { ...
    'acetylStarterH', ...
    'plastidRetainedMalonylH', ...
    'erRetainedMalonylH', ...
    'plastidKARNADPHHydrideH', ...
    'plastidENRNADHHydrideH', ...
    'plastidENRProtonH', ...
    'erKCRNADPHHydrideH', ...
    'erECRNADPHHydrideH', ...
    'erECRProtonH', ...
    'cer3AldehydeHydrideH'};

S = struct;
S.routeFractions = struct('C16',f16,'C18',f18);
S.routeDefinition = ['Fractions weight C16 versus C18 plastid-export ', ...
    'pathways before ER elongation.'];

T = struct('alkaneCarbonNumber',[],'alkaneFormula','', ...
    'alkaneHydrogenCount',[],'acylPrecursorCarbonNumber',[], ...
    'acylPrecursorLabel','', 'totalElongationCycles',[], ...
    'routeSpecific',struct,'mixed',struct);
S.compounds = repmat(T,1,numel(n));

for i = 1:numel(n)
    nc = n(i);
    alkaneH = 2*nc+2;
    precursorC = nc+1;

    R16 = routeStoichiometry(nc,16);
    R18 = routeStoichiometry(nc,18);
    M = mixRoutes(R16,R18,f16,f18,alkaneH,sourceClasses);

    C = T;
    C.alkaneCarbonNumber = nc;
    C.alkaneFormula = sprintf('C%dH%d',nc,alkaneH);
    C.alkaneHydrogenCount = alkaneH;
    C.acylPrecursorCarbonNumber = precursorC;
    C.acylPrecursorLabel = sprintf('C%d acyl precursor',precursorC);
    C.totalElongationCycles = (nc-1)/2;
    C.routeSpecific = struct('C16',R16,'C18',R18);
    C.mixed = M;
    S.compounds(i) = C;
end

% Summary table.
N = numel(S.compounds);
alkaneC=zeros(N,1); alkaneH=zeros(N,1); precursorC=zeros(N,1);
totalCycles=zeros(N,1); pCycles=zeros(N,1); erCycles=zeros(N,1);
acetyl=zeros(N,1); pMal=zeros(N,1); erMal=zeros(N,1);
pKAR=zeros(N,1); pENRH=zeros(N,1); pENRP=zeros(N,1);
erKCR=zeros(N,1); erECRH=zeros(N,1); erECRP=zeros(N,1);
cer3=zeros(N,1); totalH=zeros(N,1);

for i=1:N
    C=S.compounds(i); H=C.mixed.counts;
    alkaneC(i)=C.alkaneCarbonNumber;
    alkaneH(i)=C.alkaneHydrogenCount;
    precursorC(i)=C.acylPrecursorCarbonNumber;
    totalCycles(i)=C.totalElongationCycles;
    pCycles(i)=C.mixed.plastidElongationCycles;
    erCycles(i)=C.mixed.erElongationCycles;
    acetyl(i)=H.acetylStarterH;
    pMal(i)=H.plastidRetainedMalonylH;
    erMal(i)=H.erRetainedMalonylH;
    pKAR(i)=H.plastidKARNADPHHydrideH;
    pENRH(i)=H.plastidENRNADHHydrideH;
    pENRP(i)=H.plastidENRProtonH;
    erKCR(i)=H.erKCRNADPHHydrideH;
    erECRH(i)=H.erECRNADPHHydrideH;
    erECRP(i)=H.erECRProtonH;
    cer3(i)=H.cer3AldehydeHydrideH;
    totalH(i)=H.totalH;
end

S.summaryTable = table(alkaneC,alkaneH,precursorC,totalCycles,pCycles,erCycles, ...
    acetyl,pMal,erMal,pKAR,pENRH,pENRP,erKCR,erECRH,erECRP,cer3,totalH, ...
    'VariableNames',{'AlkaneC','AlkaneH','AcylPrecursorC', ...
    'TotalElongationCycles','PlastidElongationCycles','ERElongationCycles', ...
    'AcetylStarterH','PlastidRetainedMalonylH','ERRetainedMalonylH', ...
    'PlastidKARNADPHHydrideH','PlastidENRNADHHydrideH','PlastidENRProtonH', ...
    'ERKCRNADPHHydrideH','ERECRNADPHHydrideH','ERECRProtonH', ...
    'CER3AldehydeHydrideH','TotalH'});

S.sourceClasses = sourceClasses;
S.reactionClasses = { ...
    'plastidKAR','plastidENRHydride','plastidENRProton', ...
    'erKCR','erECRHydride','erECRProton','cer3AldehydeHydride'};
S.isotopeFractionationApplied = false;
S.isotopeCompositionApplied = false;
S.modelRole = ['reaction-explicit hydrogen-atom stoichiometry for plant ', ...
    'odd-chain n-alkane biosynthesis'];
end

function R = routeStoichiometry(alkaneC,exportC)
precursorC = alkaneC+1;
if ~ismember(exportC,[16 18])
    error('Only C16 and C18 plastid-export routes are supported.');
end
if precursorC < exportC
    error('Alkane C%d is incompatible with C%d export.',alkaneC,exportC);
end

nP=(exportC-2)/2;
nER=(precursorC-exportC)/2;
nTotal=(alkaneC-1)/2;
if nER<0 || abs(nER-round(nER))>1e-12
    error('Invalid ER elongation-cycle count.');
end

H=struct;
H.acetylStarterH=3;
H.plastidRetainedMalonylH=nP;
H.erRetainedMalonylH=nER;
H.plastidKARNADPHHydrideH=nP;
H.plastidENRNADHHydrideH=nP;
H.plastidENRProtonH=nP;
H.erKCRNADPHHydrideH=nER;
H.erECRNADPHHydrideH=nER;
H.erECRProtonH=nER;
H.cer3AldehydeHydrideH=1;
H.totalH = sumSourceCounts(H);

if abs(H.totalH-(2*alkaneC+2))>1e-12
    error('Internal H-atom bookkeeping failed for C%d.',alkaneC);
end

R=struct;
R.plastidExportCarbonNumber=exportC;
R.plastidElongationCycles=nP;
R.erElongationCycles=nER;
R.totalElongationCycles=nTotal;
R.counts=H;
R.fractions=divideCounts(H,H.totalH);
end

function M = mixRoutes(R16,R18,f16,f18,totalH,sourceClasses)
M=struct;
M.routeFractions=struct('C16',f16,'C18',f18);
M.plastidElongationCycles=f16*R16.plastidElongationCycles+f18*R18.plastidElongationCycles;
M.erElongationCycles=f16*R16.erElongationCycles+f18*R18.erElongationCycles;
H=struct;
for j=1:numel(sourceClasses)
    nm=sourceClasses{j};
    H.(nm)=f16*R16.counts.(nm)+f18*R18.counts.(nm);
end
H.totalH=sumSourceCounts(H);
if abs(H.totalH-totalH)>1e-12
    error('Route-weighted H count does not equal final alkane H count.');
end
M.counts=H;
M.fractions=divideCounts(H,totalH);
end

function total = sumSourceCounts(H)
fn=fieldnames(H);
total=0;
for i=1:numel(fn)
    if ~strcmp(fn{i},'totalH')
        total=total+H.(fn{i});
    end
end
end

function F = divideCounts(H,totalH)
F=struct;
fn=fieldnames(H);
for i=1:numel(fn)
    nm=fn{i};
    if strcmp(nm,'totalH')
        F.totalH=1;
    else
        F.(nm)=H.(nm)/totalH;
    end
end
end
