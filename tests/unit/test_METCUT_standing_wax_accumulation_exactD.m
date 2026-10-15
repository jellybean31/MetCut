function R = test_METCUT_standing_wax_accumulation_exactD()
fprintf('\n============================================================\n');
fprintf('METCUT UNIT: STANDING-WAX ACCUMULATION EXACT-D\n');
fprintf('============================================================\n');

C=[27 29 31 33]; t=0:0.25:20; nT=numel(t); details=struct;

T=wax_accumulation_default_traits();
T.productionMode="continuous"; T.lossRatePerDay=0; T.initialMass=0;
Y=repmat([-180;-160;-140;-120],1,nT);
A=wax_accumulation_dynamic_exactD(t,Y,C,T);
details.constantResidual=zeros(4,1);
for h=1:4
    idx=isfinite(A.poolDeltaD_permil(h,:));
    details.constantResidual(h)=max(abs(A.poolDeltaD_permil(h,idx)-Y(h,1)));
end
assert(max(details.constantResidual)<1e-9);

T=wax_accumulation_default_traits();
T.productionMode="developmentGated"; T.developmentWindowDays=5; T.lossRatePerDay=0;
Y=repmat(linspace(-200,-50,nT),4,1);
A=wax_accumulation_dynamic_exactD(t,Y,C,T);
idx=find(t>6);
details.lockSpread=max(A.poolDeltaD_permil(:,idx),[],2)-min(A.poolDeltaD_permil(:,idx),[],2);
assert(max(abs(details.lockSpread))<1e-10);

T=wax_accumulation_default_traits();
T.productionMode="developmentGated"; T.developmentWindowDays=0;
T.initialMass=10; T.initialDeltaD_permil=[-190;-170;-150;-130]; T.lossRatePerDay=0.1;
Y=repmat([-100;-100;-100;-100],1,nT);
A=wax_accumulation_dynamic_exactD(t,Y,C,T);
details.lossIsoResidual=max(abs(A.poolDeltaD_permil(:,end)-T.initialDeltaD_permil(:)));
expectedMass=10*exp(-0.1*t(end));
details.lossMassResidual=max(abs(A.mass(:,end)-expectedMass));
assert(details.lossIsoResidual<1e-9 && details.lossMassResidual<1e-10);

t2=[0 1];
T=wax_accumulation_default_traits();
T.productionMode="continuous"; T.peakProductionRate_perDay=1;
T.initialMass=1; T.initialDeltaD_permil=-100; T.lossRatePerDay=0;
A=wax_accumulation_dynamic_exactD(t2,zeros(4,2),C,T);
Rstd=T.Rstd_VSMOW;
xA=localD2X(-100,Rstd); xB=localD2X(0,Rstd);
expected=localX2D(0.5*(xA+xB),Rstd);
details.mixResidual=max(abs(A.poolDeltaD_permil(:,end)-expected));
assert(details.mixResidual<1e-9);

T=wax_accumulation_default_traits();
T.productionMode="developmentPlusMaintenance"; T.developmentWindowDays=2;
T.initialMass=2; T.initialDeltaD_permil=-150;
T.matureReplacementFractionPerDay=0.01; T.lossRatePerDay=0.01;
A=wax_accumulation_dynamic_exactD(t,-150*ones(4,nT),C,T);
idx=isfinite(A.poolDeltaD_permil);
details.maintenanceResidual=max(abs(A.poolDeltaD_permil(idx)+150));
assert(details.maintenanceResidual<1e-9);

R=struct('name','test_METCUT_standing_wax_accumulation_exactD','pass',true,'details',details);

fprintf('Constant-product residual:       %.3e per mil\n',max(details.constantResidual));
fprintf('Post-window lock spread:         %.3e per mil\n',max(abs(details.lockSpread)));
fprintf('Loss-only isotope residual:      %.3e per mil\n',details.lossIsoResidual);
fprintf('Equal-mass mix residual:         %.3e per mil\n',details.mixResidual);
fprintf('Maintenance residual:            %.3e per mil\n',details.maintenanceResidual);
fprintf('PASS\n');
fprintf('============================================================\n');
end

function x=localD2X(d,Rstd)
R=Rstd*(1+d/1000); x=R/(1+R);
end
function d=localX2D(x,Rstd)
R=x/(1-x); d=1000*(R/Rstd-1);
end
