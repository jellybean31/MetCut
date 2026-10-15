function W = wax_combine_hydrogen_classes_exactD(stoich,hydrogenPools)
%WAX_COMBINE_HYDROGEN_CLASSES_EXACTD
% Generic exact-D combination of H classes into whole-molecule alkane dD.
%
% This is the reaction-compatible generalization of wax_source_only_exactD.
% All upstream source-pool formation, exchange, and reaction-incorporation
% effects must already be represented in hydrogenPools.classes.
%
% Exact mixing:
%   x_wax = sum_j n_j*x_j / N_H
%
% where x is D atom fraction. Delta notation is calculated only after the
% atom-fraction mixture is complete.

if ~isstruct(stoich) || ~isfield(stoich,'compounds') || ~isfield(stoich,'sourceClasses')
    error('stoich must contain compounds and sourceClasses.');
end
if ~isstruct(hydrogenPools) || ~isfield(hydrogenPools,'classes') || ...
        ~isfield(hydrogenPools,'sourceClassNames') || ~isfield(hydrogenPools,'R_VSMOW')
    error('hydrogenPools must contain classes, sourceClassNames, and R_VSMOW.');
end

sourceNames=stoich.sourceClasses;
if ~isequal(sourceNames,hydrogenPools.sourceClassNames)
    error('Stoichiometry and H-pool class names do not match.');
end
Rstd=hydrogenPools.R_VSMOW;
if isfield(hydrogenPools,'timeDays')
    timeDays=hydrogenPools.timeDays(:).';
else
    timeDays=1:numel(hydrogenPools.classes.(sourceNames{1}).DAtomFraction);
end
nTime=numel(timeDays); nSource=numel(sourceNames); nCompound=numel(stoich.compounds);

xSource=zeros(nSource,nTime);
deltaSource=zeros(nSource,nTime);
for j=1:nSource
    nm=sourceNames{j};
    if ~isfield(hydrogenPools.classes,nm)
        error('Missing hydrogenPools.classes.%s.',nm);
    end
    C=hydrogenPools.classes.(nm);
    x=C.DAtomFraction(:).'; d=C.deltaD(:).';
    if numel(x)~=nTime || numel(d)~=nTime
        error('H class %s does not match time dimension.',nm);
    end
    xSource(j,:)=x; deltaSource(j,:)=d;
end

template=struct('alkaneCarbonNumber',[],'alkaneFormula','', ...
    'alkaneHydrogenCount',[],'routeFractions',struct,'sourceCounts',struct, ...
    'sourceFractions',struct,'sourceCountVector',[],'sourceFractionVector',[], ...
    'DAtomFraction',[],'ratio',[],'deltaD',[],'sourceDContribution',struct);
W.compounds=repmat(template,1,nCompound);
deltaDMatrix=zeros(nCompound,nTime); ratioMatrix=zeros(nCompound,nTime);
xMatrix=zeros(nCompound,nTime); sourceFractionMatrix=zeros(nCompound,nSource);

for i=1:nCompound
    C=stoich.compounds(i); H=C.mixed.counts;
    countVector=zeros(1,nSource);
    for j=1:nSource
        nm=sourceNames{j};
        if ~isfield(H,nm)
            error('Stoichiometry for C%d lacks H count %s.',C.alkaneCarbonNumber,nm);
        end
        countVector(j)=H.(nm);
    end
    NH=C.alkaneHydrogenCount;
    if abs(sum(countVector)-NH)>1e-12
        error('H counts do not sum to total H for C%d.',C.alkaneCarbonNumber);
    end
    fractionVector=countVector./NH;
    xWax=fractionVector*xSource;
    if any(~isfinite(xWax)) || any(xWax<=0) || any(xWax>=1)
        error('Calculated wax D atom fraction is nonphysical.');
    end
    RWax=xWax./(1-xWax); rWax=RWax./Rstd; dWax=1000.*(rWax-1);

    sourceCounts=struct; sourceFractions=struct; sourceDContribution=struct;
    for j=1:nSource
        nm=sourceNames{j};
        sourceCounts.(nm)=countVector(j);
        sourceFractions.(nm)=fractionVector(j);
        sourceDContribution.(nm)=fractionVector(j).*xSource(j,:);
    end
    O=template;
    O.alkaneCarbonNumber=C.alkaneCarbonNumber;
    O.alkaneFormula=C.alkaneFormula;
    O.alkaneHydrogenCount=NH;
    O.routeFractions=stoich.routeFractions;
    O.sourceCounts=sourceCounts;
    O.sourceFractions=sourceFractions;
    O.sourceCountVector=countVector;
    O.sourceFractionVector=fractionVector;
    O.DAtomFraction=xWax;
    O.ratio=rWax;
    O.deltaD=dWax;
    O.sourceDContribution=sourceDContribution;
    W.compounds(i)=O;
    deltaDMatrix(i,:)=dWax; ratioMatrix(i,:)=rWax; xMatrix(i,:)=xWax;
    sourceFractionMatrix(i,:)=fractionVector;
end

alkaneC=[W.compounds.alkaneCarbonNumber].';
tableData=array2table(sourceFractionMatrix,'VariableNames',sourceNames);
W.sourceFractionTable=addvars(tableData,alkaneC,'Before',1,'NewVariableNames','AlkaneC');
W.timeDays=timeDays;
W.R_VSMOW=Rstd;
W.sourceClassNames=sourceNames;
W.sourceDeltaDMatrix=deltaSource;
W.sourceDAtomFractionMatrix=xSource;
W.deltaDMatrix=deltaDMatrix;
W.ratioMatrix=ratioMatrix;
W.DAtomFractionMatrix=xMatrix;
W.incorporationFractionationApplied = isfield(hydrogenPools,'incorporationFractionationApplied') && ...
    hydrogenPools.incorporationFractionationApplied;
W.isotopeBookkeeping='exactDAtomFraction';
W.modelRole='whole-wax exact-D combination of resolved H classes';
end
