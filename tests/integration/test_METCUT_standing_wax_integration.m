function R = test_METCUT_standing_wax_integration()
%TEST_METCUT_STANDING_WAX_INTEGRATION
%
% METCUT integration regression for the production standing
% plant-wax accumulation operator.

thisFile=mfilename('fullpath');
integrationDir=fileparts(thisFile);
testsDir=fileparts(integrationDir);
rootDir=fileparts(testsDir);
fnDir=fullfile(rootDir,'functions');
refDir=fullfile(rootDir,'data','reference');

originalPath=path;
cleanup=onCleanup(@() path(originalPath)); 
addpath(fnDir,'-begin');
rehash

hydPath=fullfile(refDir, ...
    'hydraulic_checkpoint_76cell_8day_hourly_corrected_mass_conservative_storagefixed.mat');
assert(isfile(hydPath),'Hydraulic reference missing.');

[M,~]=run_full_model_from_hydraulic_checkpoint(hydPath);

assert(isfield(M.biochemistry,'waxStageExpanded'), ...
    'waxStageExpanded source missing.');

Wsrc=M.biochemistry.waxStageExpanded;
D=double(Wsrc.deltaDMatrix);
C=double(Wsrc.alkaneCarbonNumber(:));
if size(D,1)~=numel(C) && size(D,2)==numel(C), D=D.'; end

idx29=find(C==29,1);
assert(~isempty(idx29),'C29 missing from waxStageExpanded.');

anchorC29=D(idx29,1);
instantaneousFinalC29=D(idx29,end);

Q=load(hydPath,'timeCoarse');
timeDays=double(Q.timeCoarse(:));

Tdefault=wax_accumulation_default_traits();

W01ref_d=14;
W02support_per_d=[0.0009 0.0109];
W02ref_per_d=sqrt(prod(W02support_per_d));

schedulePath=which('wax_accumulation_production_schedule');
lossField=resolveLossFieldLocal(string(fileread(schedulePath)),Tdefault, ...
    timeDays,numel(C),W02ref_per_d);

T0=Tdefault;
T0.(char(lossField))=W02ref_per_d;

S0=wax_accumulation_production_schedule(timeDays,numel(C),T0);
peakProduction=max(double(S0.productionRate(:)));
assert(abs(peakProduction-1)<=1e-12,'Peak-production-rate regression.');

baseInitialMass=W01ref_d*peakProduction;
assert(abs(baseInitialMass-14)<=1e-12,'Standing-wax base initial-mass regression.');

massFactors=[0.5 1 2];
deltaOffsets=[-40 0 40];

finalStanding=NaN(3,3);

for im=1:3
    for io=1:3
        T=T0;
        T.initialMass=baseInitialMass*massFactors(im);
        T.initialDeltaD_permil=anchorC29+deltaOffsets(io);

        Mx=M;
        before=double(Mx.biochemistry.waxStageExpanded.deltaDMatrix);

        Mx=attach_wax_accumulation_to_model(Mx,T,'waxStageExpanded');
        A=Mx.biochemistry.waxAccumulation;

        after=double(Mx.biochemistry.waxStageExpanded.deltaDMatrix);
        product=double(A.productDeltaD_permil);

        assert(max(abs(after(:)-before(:)),[],'omitnan')<=1e-12, ...
            'Instantaneous-wax preservation regression.');
        assert(max(abs(product(:)-before(:)),[],'omitnan')<=1e-12, ...
            'Standing-wax product-input identity regression.');

        CA=double(A.alkaneCarbonNumber(:));
        ia=find(CA==29,1);
        assert(~isempty(ia),'Standing-wax C29 missing.');

        pool=double(A.poolDeltaD_permil);
        mass=double(A.mass);
        Dinv=double(A.Dinventory);

        assert(all(isfinite(pool(:))) && all(isfinite(mass(:))) && ...
            all(isfinite(Dinv(:))),'Nonfinite standing-wax output.');
        assert(all(mass(:)>=-1e-14),'Negative standing-wax mass.');

        assert(abs(mass(ia,1)-T.initialMass)<=1e-12, ...
            'Standing-wax initial-mass regression.');
        assert(abs(pool(ia,1)-T.initialDeltaD_permil)<=1e-10, ...
            'Standing-wax initial-deltaD regression.');

        assert(isfield(A,'source') && isstruct(A.source) && ...
            string(A.source.field)=="waxStageExpanded" && ...
            logical(A.source.instantaneousPreserved), ...
            'Standing-wax source metadata regression.');

        finalStanding(im,io)=pool(ia,end);
    end
end

% Accepted package-local A05d-v2 standing-wax fingerprints.
expected=[ ...
    -176.693771309, -158.151870336, -139.610093294; ...
    -189.287656296, -163.949359227, -138.611177875; ...
    -199.823656963, -168.799527831, -137.775485438];

assert(max(abs(finalStanding(:)-expected(:)))<=1e-6, ...
    'Standing-wax 3x3 final-C29 fingerprint regression.');

fullEnvelope=[min(finalStanding(:)) max(finalStanding(:))];
fullRange=diff(fullEnvelope);
central=finalStanding(2,2);
centralShift=central-instantaneousFinalC29;

assert(abs(fullRange-62.0481715251491)<=1e-6, ...
    'Standing-wax full-range regression.');
assert(abs(central-(-163.949359227336))<=1e-6, ...
    'Standing-wax central C29 regression.');
assert(abs(centralShift-(-17.6046008360977))<=1e-6, ...
    'Standing-wax central-shift regression.');

R=struct;
R.pass=true;
R.lossRateTraitField=lossField;
R.W01Reference_d=W01ref_d;
R.W02Reference_per_d=W02ref_per_d;
R.baseInitialMass=baseInitialMass;
R.finalStandingC29Grid_permil=finalStanding;
R.fullStandingC29Envelope_permil=fullEnvelope;
R.fullStandingC29Range_permil=fullRange;
R.centralFinalStandingC29_permil=central;
R.centralShiftVsInstantaneous_permil=centralShift;
end

function fld=resolveLossFieldLocal(src,T,timeDays,nH,target)
L=splitlines(src);
mask=false(numel(L),1);
for i=1:numel(L)
    if contains(lower(L(i)),"lossrate") || ...
       contains(lower(L(i)),"turnover") || contains(lower(L(i)),"loss")
        mask(max(1,i-6):min(numel(L),i+6))=true;
    end
end

tok=regexp(char(join(L(mask),newline)),'T\.([A-Za-z]\w*)','tokens');
names=strings(0,1);
for i=1:numel(tok), names(end+1,1)=string(tok{i}{1}); end %#ok<AGROW>
names=unique(names,'stable');
actual=string(fieldnames(T));
names=names(ismember(names,actual));
names=names(~ismember(names, ...
    ["initialMass","initialDeltaD_permil","Rstd_VSMOW","massTolerance"]));

ok=false(numel(names),1);
for i=1:numel(names)
    old=T.(char(names(i)));
    if isnumeric(old)&&isscalar(old)&&isfinite(old)
        Tx=T;
        Tx.(char(names(i)))=target;
        S=wax_accumulation_production_schedule(timeDays,nH,Tx);
        u=unique(double(S.lossRate(isfinite(S.lossRate))));
        ok(i)=~isempty(u)&&all(abs(u-target)<=1e-12);
    end
end

assert(nnz(ok)==1,'Could not uniquely resolve standing-wax loss-rate trait.');
fld=names(ok);
end
