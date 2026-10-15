function R = test_METCUT_intracellular_water_exactD()
fprintf('\n============================================================\n');
fprintf('METCUT UNIT: INTRACELLULAR WATER EXACT-D\n');
fprintf('============================================================\n');

Rstd=155.76e-6;
ratioToX=@(r)(Rstd.*r)./(1+Rstd.*r);
xToRatio=@(x)(x./(1-x))./Rstd;
ratioToDelta=@(r)1000.*(r-1);

cfg=struct;
cfg.includePlastid=true;
cfg.includeER=false;
cfg.waterFractionOfSymplast=struct('cytosol',0.05,'plastid',0.003);
cfg.exchangeTimescale_s=struct('symplastCytosol',0.5,'cytosolPlastid',0.1);
cfg.sourcePool='symplastBulk';
cfg.sourceInterpolation='linear';
cfg.initialState='sourceEqual';
cfg.odeRelTol=1e-10;
cfg.odeAbsTol=1e-14;
cfg.odeMaxStepFraction=1;

details=struct;

t1=(0:1/48:3).'; d0=-75;
WS1=6.5+0.7.*sin(2*pi.*t1./1.2)+0.25.*cos(2*pi.*t1./0.55);
leaf1=struct('deltaDSymplastBulk',d0.*ones(size(t1)), ...
    'symplasticWaterAmount_mol_m2',WS1);
I1=intracellular_water_dynamic_exactD(t1,leaf1,cfg);

details.uniformDeviation=max([ ...
    max(abs(I1.deltaD.cytosol-d0),[],'omitnan'), ...
    max(abs(I1.deltaD.plastid-d0),[],'omitnan')]);
details.uniformWaterResidual=max([ ...
    max(abs(I1.compartmentWaterClosureResidual_mol_m2.cytosol),[],'omitnan'), ...
    max(abs(I1.compartmentWaterClosureResidual_mol_m2.plastid),[],'omitnan')]);
details.uniformDResidual=max([ ...
    max(abs(I1.compartmentDClosureResidual_molD_m2.cytosol),[],'omitnan'), ...
    max(abs(I1.compartmentDClosureResidual_molD_m2.plastid),[],'omitnan'), ...
    max(abs(I1.wholeSystemDClosureResidual_molD_m2),[],'omitnan')]);

assert(details.uniformDeviation<1e-7);
assert(details.uniformWaterResidual<1e-12);
assert(details.uniformDResidual<1e-12);

t2=(0:1/288:0.25).'; WS2=6.0;
leaf2=struct('deltaDSymplastBulk',-75.*ones(size(t2)), ...
    'symplasticWaterAmount_mol_m2',WS2.*ones(size(t2)));
cfg2=cfg;
cfg2.exchangeTimescale_s.symplastCytosol=inf;
cfg2.exchangeTimescale_s.cytosolPlastid=15;
cfg2.initialState='specified';
cfg2.initialDeltaD_permil=struct('cytosol',-50,'plastid',-150);
cfg2.odeRelTol=1e-12; cfg2.odeAbsTol=1e-18; cfg2.odeMaxStepFraction=0.25;
I2=intracellular_water_dynamic_exactD(t2,leaf2,cfg2);

WCy=cfg2.waterFractionOfSymplast.cytosol.*WS2;
WP=cfg2.waterFractionOfSymplast.plastid.*WS2;
QCyP=WP./cfg2.exchangeTimescale_s.cytosolPlastid;
xCy0=ratioToX(1-50/1000); xP0=ratioToX(1-150/1000);
Wtot=WCy+WP; xEq=(WCy.*xCy0+WP.*xP0)./Wtot;
lambda=QCyP.*(1./WCy+1./WP);
difference=(xCy0-xP0).*exp(-lambda.*t2.*86400);
xCyExpected=xEq+(WP./Wtot).*difference;
xPExpected=xEq-(WCy./Wtot).*difference;

details.relaxationResidual=max([ ...
    max(abs(I2.deltaD.cytosol(:)-ratioToDelta(xToRatio(xCyExpected))),[],'omitnan'), ...
    max(abs(I2.deltaD.plastid(:)-ratioToDelta(xToRatio(xPExpected))),[],'omitnan')]);
details.closedDResidual=max(abs(I2.totalDInventory_molD_m2-I2.totalDInventory_molD_m2(1)),[],'omitnan');
assert(details.relaxationResidual<1e-6);
assert(details.closedDResidual<1e-12);

t3=(0:1/48:3).';
source3=-70+20.*sin(2*pi.*t3./0.85)-6.*cos(2*pi.*t3./0.37);
WS3=6.2+0.9.*sin(2*pi.*t3./1.05);
leaf3=struct('deltaDSymplastBulk',source3,'symplasticWaterAmount_mol_m2',WS3);
I3=intracellular_water_dynamic_exactD(t3,leaf3,cfg);

inventoryChange=I3.totalDInventory_molD_m2(end)-I3.totalDInventory_molD_m2(1);
externalD=I3.cumulativeDTransfer_molD_m2.symplastCytosol(end);
details.variableBoundaryResidual=abs(inventoryChange-externalD);
details.variableDResidual=max([ ...
    max(abs(I3.wholeSystemDClosureResidual_molD_m2),[],'omitnan'), ...
    max(abs(I3.compartmentDClosureResidual_molD_m2.cytosol),[],'omitnan'), ...
    max(abs(I3.compartmentDClosureResidual_molD_m2.plastid),[],'omitnan')]);
details.variableWaterResidual=max([ ...
    max(abs(I3.compartmentWaterClosureResidual_mol_m2.cytosol),[],'omitnan'), ...
    max(abs(I3.compartmentWaterClosureResidual_mol_m2.plastid),[],'omitnan')]);

assert(details.variableBoundaryResidual<1e-12);
assert(details.variableDResidual<1e-12);
assert(details.variableWaterResidual<1e-12);

R=struct('name','test_METCUT_intracellular_water_exactD','pass',true,'details',details);

fprintf('Uniform isotope deviation:       %.3e per mil\n',details.uniformDeviation);
fprintf('Closed-pair analytical residual: %.3e per mil\n',details.relaxationResidual);
fprintf('Variable-forcing D residual:     %.3e mol D m^-2\n',details.variableDResidual);
fprintf('PASS\n');
fprintf('============================================================\n');
end
