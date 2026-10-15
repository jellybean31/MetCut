
function R = test_wax_accumulation_dynamic_exactD()
%TEST_WAX_ACCUMULATION_DYNAMIC_EXACTD

fprintf('\n');
fprintf('============================================================\n');
fprintf('WAX ACCUMULATION EXACT-D REGRESSION TESTS\n');
fprintf('============================================================\n');

C=[27 29 31 33];
t=0:0.25:20;
nT=numel(t);

pass=true;
details=struct;

%% 1. Constant product isotope composition must be preserved.
T=wax_accumulation_default_traits();
T.productionMode="continuous";
T.lossRatePerDay=0;
T.initialMass=0;

Y=repmat([-180;-160;-140;-120],1,nT);

A=wax_accumulation_dynamic_exactD(t,Y,C,T);

for h=1:4
    idx=isfinite(A.poolDeltaD_permil(h,:));
    err=max(abs(A.poolDeltaD_permil(h,idx)-Y(h,1)));
    details.constantProductResidual(h,1)=err;
end

details.constantProductPass=max(details.constantProductResidual)<1e-9;
pass=pass && details.constantProductPass;

%% 2. Development-gated pool must lock after production stops if loss=0.
T=wax_accumulation_default_traits();
T.productionMode="developmentGated";
T.developmentWindowDays=5;
T.lossRatePerDay=0;

Y=repmat(linspace(-200,-50,nT),4,1);

A=wax_accumulation_dynamic_exactD(t,Y,C,T);

idx=find(t>6);
lockSpread=max(A.poolDeltaD_permil(:,idx),[],2)- ...
    min(A.poolDeltaD_permil(:,idx),[],2);

details.developmentLockSpread=lockSpread;
details.developmentLockPass=max(abs(lockSpread))<1e-10;
pass=pass && details.developmentLockPass;

%% 3. Loss without production changes mass but not isotope composition.
T=wax_accumulation_default_traits();
T.productionMode="developmentGated";
T.developmentWindowDays=0;
T.initialMass=10;
T.initialDeltaD_permil=[-190;-170;-150;-130];
T.lossRatePerDay=0.1;

Y=repmat([-100;-100;-100;-100],1,nT);

A=wax_accumulation_dynamic_exactD(t,Y,C,T);

isoErr=max(abs(A.poolDeltaD_permil(:,end)- ...
    T.initialDeltaD_permil(:)));

expectedMass=10*exp(-0.1*t(end));
massErr=max(abs(A.mass(:,end)-expectedMass));

details.lossOnlyIsotopeResidual=isoErr;
details.lossOnlyMassResidual=massErr;
details.lossOnlyPass=isoErr<1e-9 && massErr<1e-10;
pass=pass && details.lossOnlyPass;

%% 4. Exact-D equal-mass mixture.
t2=[0 1];
T=wax_accumulation_default_traits();
T.productionMode="continuous";
T.peakProductionRate_perDay=1;
T.initialMass=1;
T.initialDeltaD_permil=-100;
T.lossRatePerDay=0;

Y2=zeros(4,2);

A=wax_accumulation_dynamic_exactD(t2,Y2,C,T);

Rstd=T.Rstd_VSMOW;
xA=d2x(-100,Rstd);
xB=d2x(0,Rstd);
expected=x2d(0.5*(xA+xB),Rstd);

mixErr=max(abs(A.poolDeltaD_permil(:,end)-expected));

details.equalMassExactDMixtureResidual=mixErr;
details.equalMassExactDMixturePass=mixErr<1e-9;
pass=pass && details.equalMassExactDMixturePass;

%% 5. Mature maintenance with product equal to existing pool preserves isotope.
T=wax_accumulation_default_traits();
T.productionMode="developmentPlusMaintenance";
T.developmentWindowDays=2;
T.initialMass=2;
T.initialDeltaD_permil=-150;
T.matureReplacementFractionPerDay=0.01;
T.lossRatePerDay=0.01;

Y3=-150*ones(4,nT);

A=wax_accumulation_dynamic_exactD(t,Y3,C,T);

idx=isfinite(A.poolDeltaD_permil);
maintErr=max(abs(A.poolDeltaD_permil(idx)+150));

details.maintenanceIsotopeResidual=maintErr;
details.maintenancePass=maintErr<1e-9;
pass=pass && details.maintenancePass;

R=struct;
R.details=details;
R.overall=pass;

fprintf('Constant product exact-D preservation: %d | max err %.3e per mil\n', ...
    details.constantProductPass,max(details.constantProductResidual));
fprintf('Development-gated locking: %d | max spread %.3e per mil\n', ...
    details.developmentLockPass,max(details.developmentLockSpread));
fprintf('Loss-only isotope/mass analytic test: %d | %.3e / %.3e\n', ...
    details.lossOnlyPass, ...
    details.lossOnlyIsotopeResidual, ...
    details.lossOnlyMassResidual);
fprintf('Equal-mass exact-D mixture: %d | err %.3e per mil\n', ...
    details.equalMassExactDMixturePass, ...
    details.equalMassExactDMixtureResidual);
fprintf('Mature maintenance isotope preservation: %d | err %.3e per mil\n', ...
    details.maintenancePass, ...
    details.maintenanceIsotopeResidual);
fprintf('OVERALL PASS: %d\n',R.overall);
fprintf('============================================================\n');

end

function x=d2x(d,Rstd)
R=Rstd*(1+d/1000);
x=R/(1+R);
end

function d=x2d(x,Rstd)
R=x/(1-x);
d=1000*(R/Rstd-1);
end
