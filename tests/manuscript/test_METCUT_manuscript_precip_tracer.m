function R = test_METCUT_manuscript_precip_tracer()
%TEST_METCUT_MANUSCRIPT_PRECIP_TRACER
%
% Regression for the paired precipitation
% isotope tracer experiment.
%
% Experiment:
%   baseline precipitation deltaD = -50 permil
%   pulse intervals              = 7:10
%   pulse window                 = 0.250000 -> 0.416667 d
%   positive pulse               = +100 permil anomaly
%   negative pulse               = -100 permil anomaly
%
% Hydraulics are reused.

thisFile=mfilename('fullpath');
manuscriptDir=fileparts(thisFile);
testsDir=fileparts(manuscriptDir);
rootDir=fileparts(testsDir);
fnDir=fullfile(rootDir,'functions');
refDir=fullfile(rootDir,'data','reference');

originalPath=path;
cleanup=onCleanup(@() path(originalPath)); 
addpath(fnDir,'-begin');
rehash

hydPath=fullfile(refDir, ...
    'hydraulic_checkpoint_76cell_8day_hourly_corrected_mass_conservative_storagefixed.mat');
assert(isfile(hydPath),'Package-local hydraulic reference missing.');

[M0,B0]=run_full_model_from_hydraulic_checkpoint(hydPath);

t=double(M0.timeDays(:));
nT=numel(t);
nInt=nT-1;

d0=expandIntervalLocal( ...
    M0.config.forcing.precipDeltaD_permil,nInt, ...
    'baseline precipitation deltaD');

assert(max(abs(d0+50))<=1e-12, ...
    'Tracer baseline precipitation must be -50 permil.');

pulseIdx=(7:10).';
assert(nInt>=pulseIdx(end),'Tracer pulse indices exceed run length.');

pulseStart=t(pulseIdx(1));
pulseEnd=t(pulseIdx(end)+1);
pulseMid=0.5.*(t(1:end-1)+t(2:end));
pulseCentroid=mean(pulseMid(pulseIdx));

assert(abs(pulseStart-0.25)<=1e-12, ...
    'Tracer pulse-start regression.');
assert(abs(pulseEnd-(10/24))<=1e-12, ...
    'Tracer pulse-end regression.');
assert(abs(pulseCentroid-(1/3))<=1e-12, ...
    'Tracer pulse-centroid regression.');

% All four accepted pulse intervals must carry infiltration.
if isfield(M0.hydraulics,'soil') && ...
        isfield(M0.hydraulics.soil,'JinfiltrationActualInterval_mol_m2_s')
    Jin=double(M0.hydraulics.soil.JinfiltrationActualInterval_mol_m2_s(:));
    assert(all(Jin(pulseIdx)>0), ...
        'One or more tracer pulse intervals has no infiltration.');
end

dPlus=d0;
dMinus=d0;
dPlus(pulseIdx)=dPlus(pulseIdx)+100;
dMinus(pulseIdx)=dMinus(pulseIdx)-100;

[Mplus,Bplus]=rerun_full_model_isotope_forcing(M0,dPlus);
[Mminus,Bminus]=rerun_full_model_isotope_forcing(M0,dMinus);

% Replays must not modify the hydraulic trajectory.
assert(isequaln(Mplus.hydraulics,M0.hydraulics), ...
    'Positive tracer replay modified hydraulics.');
assert(isequaln(Mminus.hydraulics,M0.hydraulics), ...
    'Negative tracer replay modified hydraulics.');

F=struct;
F.baselineDeltaDInterval_permil=d0;
F.plusDeltaDInterval_permil=dPlus;
F.minusDeltaDInterval_permil=dMinus;
F.pulseIntervalIndex=pulseIdx;

T=analyze_full_model_precip_tracer_experiment( ...
    M0,B0,Mplus,Bplus,Mminus,Bminus,F,struct('makeFigures',false));

assert(abs(T.forcing.inputOddAmplitude_permil-100)<=1e-12, ...
    'Tracer odd-input amplitude regression.');
assert(T.forcing.inputEvenMax_permil<=1e-12, ...
    'Tracer forcing is not an exact +/- mirror pair.');

rootPeak=T.response.rootSource.metrics.peakOddResponse_permil;
xylemPeak=T.response.xylemDelivered.metrics.peakOddResponse_permil;
leafPeak=T.response.leafSymplastBulk.metrics.peakOddResponse_permil;
c29Peak=T.wax.response.waxC29.metrics.peakOddResponse_permil;
c29Tail=T.wax.response.waxC29.metrics.tailAreaFractionAfterPulse;

rootPeakTime=T.response.rootSource.metrics.peakTime_days;
xylemPeakTime=T.response.xylemDelivered.metrics.peakTime_days;
leafPeakTime=T.response.leafSymplastBulk.metrics.peakTime_days;
c29PeakTime=T.wax.response.waxC29.metrics.peakTime_days;

assert(abs(rootPeak-14.5844)<=5e-4, ...
    'Root-source tracer-peak regression.');
assert(abs(xylemPeak-10.2451571232)<=2e-6, ...
    'Xylem tracer-peak regression.');
assert(abs(leafPeak-2.9405)<=5e-4, ...
    'Bulk-leaf tracer-peak regression.');
assert(abs(c29Peak-2.57612)<=8e-5, ...
    'C29 tracer-peak regression.');
assert(abs(c29Tail-0.99474)<=8e-5, ...
    'C29 finite-window tail-fraction regression.');

assert(abs(rootPeakTime-7.25)<=1/24+1e-12, ...
    'Root-source tracer peak-time regression.');
assert(abs(xylemPeakTime-5.79166666666667)<=1/24+1e-12, ...
    'Xylem tracer peak-time regression.');
assert(abs(leafPeakTime-6.29166666666667)<=1/24+1e-12, ...
    'Bulk-leaf tracer peak-time regression.');
assert(abs(c29PeakTime-6.29166666666667)<=1/24+1e-12, ...
    'C29 tracer peak-time regression.');

% Multihomologue response must remain tightly clustered and ordered.
C=double(T.wax.homologSummary.AlkaneC(:));
P=double(T.wax.homologSummary.PeakOddResponse_permil(:));
[Cs,ord]=sort(C);
Ps=P(ord);

assert(isequal(Cs(:).',[27 29 31 33]), ...
    'Tracer homolog set regression.');
assert(all(diff(Ps)>0), ...
    'Tracer homolog-peak ordering regression.');

expectedHomologPeak=[2.5755 2.57612 2.5767 2.5772];
assert(max(abs(Ps(:).'-expectedHomologPeak))<=1e-3, ...
    'Tracer multihomologue response fingerprint regression.');

% Exact-D closure must remain within accepted integration tolerances.
checkClosureLocal(B0,'baseline');
checkClosureLocal(Bplus,'positive tracer');
checkClosureLocal(Bminus,'negative tracer');

% Gross response must be attenuated relative to the 100-permil input.
assert(rootPeak>0 && rootPeak<100, ...
    'Root tracer response outside expected physical range.');
assert(xylemPeak>0 && xylemPeak<rootPeak, ...
    'Xylem tracer attenuation regression.');
assert(leafPeak>0 && leafPeak<xylemPeak, ...
    'Leaf tracer attenuation regression.');
assert(c29Peak>0 && c29Peak<leafPeak, ...
    'Plant-wax tracer attenuation regression.');

fprintf('\n');
fprintf('============================================================\n');
fprintf('METCUT PRECIPITATION-TRACER REGRESSION\n');
fprintf('============================================================\n');
fprintf('Pulse intervals:                 %s\n',mat2str(pulseIdx.'));
fprintf('Pulse start/end:                 %.9f -> %.9f d\n',pulseStart,pulseEnd);
fprintf('Root-source peak:                %.12g permil\n',rootPeak);
fprintf('Xylem-delivered peak:            %.12g permil\n',xylemPeak);
fprintf('Bulk-symplast peak:              %.12g permil\n',leafPeak);
fprintf('Instantaneous C29 peak:          %.12g permil\n',c29Peak);
fprintf('C29 tail fraction after pulse:   %.12g\n',c29Tail);
fprintf('Homolog peaks C27/C29/C31/C33:  %s permil\n',mat2str(Ps.',12));
fprintf('PASS\n');
fprintf('============================================================\n');

R=struct;
R.pass=true;
R.pulseIntervalIndex=pulseIdx;
R.pulseStart_days=pulseStart;
R.pulseEnd_days=pulseEnd;
R.pulseCentroid_days=pulseCentroid;
R.rootSourcePeak_permil=rootPeak;
R.xylemDeliveredPeak_permil=xylemPeak;
R.leafSymplastBulkPeak_permil=leafPeak;
R.C29Peak_permil=c29Peak;
R.C29TailFraction=c29Tail;
R.homologCarbonNumber=Cs;
R.homologPeak_permil=Ps;
R.analysis=T;
end

% =========================================================================
function x=expandIntervalLocal(x,nInt,label)
if isscalar(x)
    x=repmat(double(x),nInt,1);
else
    x=double(x(:));
end
assert(numel(x)==nInt && all(isfinite(x)), ...
    '%s must resolve to nTime-1 finite values.',label);
end

% =========================================================================
function checkClosureLocal(B,label)
assert(isfield(B,'closure'),'%s budget lacks closure diagnostics.',label);

assert(B.closure.maxSoilDResidual_molD_m2<1e-10, ...
    '%s soil exact-D closure regression.',label);
assert(B.closure.maxXylemDResidual_molD_m2<1e-10, ...
    '%s xylem exact-D closure regression.',label);
assert(B.closure.maxLeafDResidual_molD_m2<1e-9, ...
    '%s leaf exact-D closure regression.',label);
assert(B.closure.maxIntracellularDResidual_molD_m2<1e-12, ...
    '%s intracellular exact-D closure regression.',label);
end
