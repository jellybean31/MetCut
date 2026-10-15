function T = analyze_full_model_precip_tracer_experiment( ...
    Mbase,Bbase,Mplus,Bplus,Mminus,Bminus,F,cfg)
%ANALYZE_FULL_MODEL_PRECIP_TRACER_EXPERIMENT
% Differential positive/negative precipitation-isotope tracer analysis.
%
% Hydraulics are identical in all three runs. The analysis therefore uses
% tracer-minus-baseline anomalies to isolate propagation of the imposed
% precipitation isotope signal from meteorological/physiological variability
% that is common to all runs.
%
% Positive case:
%   baseline -50 per mil -> +50 per mil pulse -> baseline
%
% Negative mirror case:
%   baseline -50 per mil -> -150 per mil pulse -> baseline
%
% The two pulse amplitudes are therefore +/-100 per mil.
%
% Differential decomposition
% --------------------------
% Let:
%
%   p(t) = delta_plus(t)  - delta_baseline(t)
%   m(t) = delta_minus(t) - delta_baseline(t)
%
% Then:
%
%   odd response  = 0.5 * [p(t) - m(t)]
%   even response = 0.5 * [p(t) + m(t)]
%
% Perfect mirror symmetry gives even response = 0.
%
% Because deltaD = 1000*(R/Rstd - 1), deltaD anomalies are exactly
% proportional to isotope-ratio anomalies. We therefore report anomaly
% amplitudes in per mil while retaining exact-D forward calculations.

if nargin<8 || isempty(cfg), cfg=struct; end
cfg=setDefaultsLocal(cfg);

t=Mbase.timeDays(:);
nT=numel(t);
nInt=nT-1;

validateRunDimensionsLocal(Mbase,Mplus,Mminus,nT);

d0=F.baselineDeltaDInterval_permil(:);
dp=F.plusDeltaDInterval_permil(:);
dm=F.minusDeltaDInterval_permil(:);

if any([numel(d0),numel(dp),numel(dm)]~=nInt)
    error('Tracer forcing vectors must each contain nTime-1 values.');
end

pulseIdx=F.pulseIntervalIndex(:);
if isempty(pulseIdx)
    error('At least one pulse interval is required.');
end

dt=diff(t);
tmid=0.5.*(t(1:end-1)+t(2:end));

inputPlus=dp-d0;
inputMinus=dm-d0;
inputOdd=0.5.*(inputPlus-inputMinus);
inputEven=0.5.*(inputPlus+inputMinus);

inputAmplitude=max(inputOdd);
if ~isfinite(inputAmplitude) || inputAmplitude<=0
    error('Positive odd tracer amplitude could not be determined.');
end

inputArea=sum(max(inputOdd,0).*dt);
if inputArea<=0
    error('Tracer pulse has zero positive time-integrated input.');
end

inputCentroid=sum(tmid.*max(inputOdd,0).*dt)./inputArea;
pulseStart=min(t(pulseIdx));
pulseEnd=max(t(pulseIdx+1));

%% ========================================================================
% Primary state responses
% ========================================================================

stateNames={ ...
    'soilBulk', ...
    'rootSource', ...
    'xylemDelivered', ...
    'leafApoplast', ...
    'leafSymplastBulk', ...
    'evaporationSite', ...
    'cytosol', ...
    'plastid', ...
    'preReactionC29Mix', ...
    'waxC29'};

nS=numel(stateNames);
response=struct;
rows=cell(nS,22);

for i=1:nS
    nm=stateNames{i};

    b=Bbase.states.(nm).deltaD_permil(:);
    p=Bplus.states.(nm).deltaD_permil(:)-b;
    m=Bminus.states.(nm).deltaD_permil(:)-b;

    odd=0.5.*(p-m);
    even=0.5.*(p+m);

    R=responseMetricsLocal( ...
        t,p,m,odd,even,inputAmplitude,inputArea, ...
        pulseStart,pulseEnd,inputCentroid,cfg);

    response.(nm)=struct;
    response.(nm).plusAnomaly_permil=p;
    response.(nm).minusAnomaly_permil=m;
    response.(nm).oddResponse_permil=odd;
    response.(nm).evenResponse_permil=even;
    response.(nm).normalizedOddResponse=odd./inputAmplitude;
    response.(nm).metrics=R;

    rows(i,:)={ ...
        nm, ...
        R.peakPlusOriented_permil, ...
        R.peakMinusOriented_permil, ...
        R.peakOddResponse_permil, ...
        R.peakTransferFractionOfRain, ...
        R.peakTime_days, ...
        R.peakLagFromPulseCentroid_days, ...
        R.onsetTime05_days, ...
        R.onsetLagFromPulseStart_days, ...
        R.centroidTime_days, ...
        R.centroidLagFromPulse_days, ...
        R.positiveArea_permil_days, ...
        R.netArea_permil_days, ...
        R.areaGainVsRain, ...
        R.netAreaGainVsRain, ...
        R.tailAreaFractionAfterPulse, ...
        R.finalResidualFractionOfPeak, ...
        R.negativeLobeAreaFraction, ...
        R.symmetryMaxResidualFractionOfInput, ...
        R.symmetryRMSResidualFractionOfInput, ...
        R.plusVsOrientedMinusCorrelation, ...
        R.windowLikelyTruncated};
end

stateResponseSummary=cell2table(rows,'VariableNames',{ ...
    'State', ...
    'PeakPlusOriented_permil', ...
    'PeakMinusOriented_permil', ...
    'PeakOddResponse_permil', ...
    'PeakTransferFractionOfRain', ...
    'PeakTime_days', ...
    'PeakLagFromPulseCentroid_days', ...
    'OnsetTime05_days', ...
    'OnsetLagFromPulseStart_days', ...
    'CentroidTime_days', ...
    'CentroidLagFromPulse_days', ...
    'PositiveArea_permil_days', ...
    'NetArea_permil_days', ...
    'StateResponseAreaGainVsRain', ...
    'NetStateResponseAreaGainVsRain', ...
    'TailAreaFractionAfterPulse', ...
    'FinalResidualFractionOfPeak', ...
    'NegativeLobeAreaFraction', ...
    'SymmetryMaxResidualFractionOfInput', ...
    'SymmetryRMSResidualFractionOfInput', ...
    'PlusVsOrientedMinusCorrelation', ...
    'WindowLikelyTruncated'});

%% ========================================================================
% Causal transfer chain to wax
% ========================================================================

causalStates={ ...
    'precipitationPulse', ...
    'soilBulk', ...
    'rootSource', ...
    'xylemDelivered', ...
    'leafApoplast', ...
    'leafSymplastBulk', ...
    'cytosol', ...
    'plastid', ...
    'preReactionC29Mix', ...
    'waxC29'};

nC=numel(causalStates);

peak=zeros(nC,1);
area=zeros(nC,1);
centroid=zeros(nC,1);

peak(1)=inputAmplitude;
area(1)=inputArea;
centroid(1)=inputCentroid;

for i=2:nC
    R=response.(causalStates{i}).metrics;
    peak(i)=R.peakOddResponse_permil;
    area(i)=R.positiveArea_permil_days;
    centroid(i)=R.centroidTime_days;
end

cRows=cell(nC-1,10);

for i=1:nC-1
    if abs(peak(i))>0
        peakGain=peak(i+1)/peak(i);
    else
        peakGain=NaN;
    end

    if abs(area(i))>0
        areaGain=area(i+1)/area(i);
    else
        areaGain=NaN;
    end

    if isfinite(centroid(i)) && isfinite(centroid(i+1))
        incrementalLag=centroid(i+1)-centroid(i);
    else
        incrementalLag=NaN;
    end

    Rout=response.(causalStates{i+1}).metrics;

    cRows(i,:)={ ...
        causalStates{i},causalStates{i+1}, ...
        peakGain,areaGain,incrementalLag, ...
        Rout.peakTransferFractionOfRain, ...
        Rout.centroidLagFromPulse_days, ...
        Rout.tailAreaFractionAfterPulse, ...
        Rout.symmetryMaxResidualFractionOfInput, ...
        Rout.windowLikelyTruncated};
end

causalTransferSummary=cell2table(cRows,'VariableNames',{ ...
    'FromStage','ToStage', ...
    'PeakTracerGain_OutOverIn', ...
    'StateResponseAreaGain_OutOverIn', ...
    'IncrementalCentroidLag_days', ...
    'OutputPeakTransferFractionOfRain', ...
    'OutputCentroidLagFromRain_days', ...
    'OutputTailAreaFractionAfterPulse', ...
    'OutputSymmetryMaxResidualFractionOfInput', ...
    'OutputWindowLikelyTruncated'});

%% ========================================================================
% Resolved soil-cell tracer propagation
% ========================================================================

dbase=Mbase.soilIsotopes.deltaD_permil;
dplus=Mplus.soilIsotopes.deltaD_permil-dbase;
dminus=Mminus.soilIsotopes.deltaD_permil-dbase;

oddSoil=0.5.*(dplus-dminus);
evenSoil=0.5.*(dplus+dminus);

[nCell,nTs]=size(oddSoil);
if nTs~=nT
    error('Resolved soil isotope trajectory has wrong time dimension.');
end

z=Mbase.grid.zCenter_m(:);
zone=Mbase.grid.zoneIndex(:);

if numel(z)~=nCell || numel(zone)~=nCell
    error('Grid metadata dimensions do not match soil isotope cells.');
end

soilRows=cell(nCell,10);

for c=1:nCell
    p=dplus(c,:).';
    m=dminus(c,:).';
    odd=oddSoil(c,:).';
    even=evenSoil(c,:).';

    R=responseMetricsLocal( ...
        t,p,m,odd,even,inputAmplitude,inputArea, ...
        pulseStart,pulseEnd,inputCentroid,cfg);

    soilRows(c,:)={ ...
        c,z(c),zone(c), ...
        R.peakOddResponse_permil, ...
        R.peakTransferFractionOfRain, ...
        R.peakTime_days, ...
        R.centroidLagFromPulse_days, ...
        R.tailAreaFractionAfterPulse, ...
        R.symmetryMaxResidualFractionOfInput, ...
        R.windowLikelyTruncated};
end

soilCellSummary=cell2table(soilRows,'VariableNames',{ ...
    'CellIndex','DepthCenter_m','RootingZone', ...
    'PeakOddResponse_permil','PeakTransferFractionOfRain', ...
    'PeakTime_days','CentroidLagFromPulse_days', ...
    'TailAreaFractionAfterPulse', ...
    'SymmetryMaxResidualFractionOfInput', ...
    'WindowLikelyTruncated'});

%% ========================================================================
% Soil-zone exact inventory-weighted tracer response
% ========================================================================

zones=unique(zone(:).');
zoneRows=cell(numel(zones),10);
soilZoneResponse=struct;

for i=1:numel(zones)
    iz=zones(i);

    dB=zoneBulkDeltaLocal(Mbase,iz);
    dPz=zoneBulkDeltaLocal(Mplus,iz)-dB;
    dMz=zoneBulkDeltaLocal(Mminus,iz)-dB;

    odd=0.5.*(dPz-dMz);
    even=0.5.*(dPz+dMz);

    R=responseMetricsLocal( ...
        t,dPz,dMz,odd,even,inputAmplitude,inputArea, ...
        pulseStart,pulseEnd,inputCentroid,cfg);

    fn=sprintf('zone%d',iz);
    soilZoneResponse.(fn).oddResponse_permil=odd;
    soilZoneResponse.(fn).metrics=R;

    zoneRows(i,:)={ ...
        iz,sum(zone==iz), ...
        min(z(zone==iz)),max(z(zone==iz)), ...
        R.peakOddResponse_permil, ...
        R.peakTransferFractionOfRain, ...
        R.peakTime_days, ...
        R.centroidLagFromPulse_days, ...
        R.tailAreaFractionAfterPulse, ...
        R.windowLikelyTruncated};
end

soilZoneSummary=cell2table(zoneRows,'VariableNames',{ ...
    'RootingZone','Ncells','MinDepthCenter_m','MaxDepthCenter_m', ...
    'PeakOddResponse_permil','PeakTransferFractionOfRain', ...
    'PeakTime_days','CentroidLagFromPulse_days', ...
    'TailAreaFractionAfterPulse','WindowLikelyTruncated'});

%% ========================================================================
% Wax homolog tracer response
% ========================================================================

chains=Mbase.biochemistry.alkaneCarbonNumber(:);
waxRows=cell(numel(chains),12);
waxResponse=struct;

for i=1:numel(chains)
    nm=sprintf('waxC%d',chains(i));

    b=Bbase.states.(nm).deltaD_permil(:);
    p=Bplus.states.(nm).deltaD_permil(:)-b;
    m=Bminus.states.(nm).deltaD_permil(:)-b;
    odd=0.5.*(p-m);
    even=0.5.*(p+m);

    R=responseMetricsLocal( ...
        t,p,m,odd,even,inputAmplitude,inputArea, ...
        pulseStart,pulseEnd,inputCentroid,cfg);

    waxResponse.(nm).oddResponse_permil=odd;
    waxResponse.(nm).evenResponse_permil=even;
    waxResponse.(nm).metrics=R;

    waxRows(i,:)={ ...
        chains(i),R.peakOddResponse_permil,R.peakTransferFractionOfRain, ...
        R.peakTime_days,R.centroidLagFromPulse_days, ...
        R.areaGainVsRain,R.tailAreaFractionAfterPulse, ...
        R.finalResidualFractionOfPeak, ...
        R.symmetryMaxResidualFractionOfInput, ...
        R.symmetryRMSResidualFractionOfInput, ...
        R.plusVsOrientedMinusCorrelation, ...
        R.windowLikelyTruncated};
end

waxHomologSummary=cell2table(waxRows,'VariableNames',{ ...
    'AlkaneC','PeakOddResponse_permil','PeakTransferFractionOfRain', ...
    'PeakTime_days','CentroidLagFromPulse_days', ...
    'StateResponseAreaGainVsRain','TailAreaFractionAfterPulse', ...
    'FinalResidualFractionOfPeak', ...
    'SymmetryMaxResidualFractionOfInput', ...
    'SymmetryRMSResidualFractionOfInput', ...
    'PlusVsOrientedMinusCorrelation', ...
    'WindowLikelyTruncated'});

%% ========================================================================
% Closure summary
% ========================================================================

closureRows={ ...
    'baseline', ...
        Bbase.closure.maxSoilDResidual_molD_m2, ...
        Bbase.closure.maxXylemDResidual_molD_m2, ...
        Bbase.closure.maxLeafDResidual_molD_m2, ...
        Bbase.closure.maxIntracellularDResidual_molD_m2; ...
    'positivePulse', ...
        Bplus.closure.maxSoilDResidual_molD_m2, ...
        Bplus.closure.maxXylemDResidual_molD_m2, ...
        Bplus.closure.maxLeafDResidual_molD_m2, ...
        Bplus.closure.maxIntracellularDResidual_molD_m2; ...
    'negativePulse', ...
        Bminus.closure.maxSoilDResidual_molD_m2, ...
        Bminus.closure.maxXylemDResidual_molD_m2, ...
        Bminus.closure.maxLeafDResidual_molD_m2, ...
        Bminus.closure.maxIntracellularDResidual_molD_m2};

closureSummary=cell2table(closureRows,'VariableNames',{ ...
    'Run','MaxSoilDResidual_molD_m2','MaxXylemDResidual_molD_m2', ...
    'MaxLeafDResidual_molD_m2','MaxIntracellularDResidual_molD_m2'});

%% ========================================================================
% Output
% ========================================================================

T=struct;
T.version='2026-08-19_precipitation_isotope_differential_tracer_v1';
T.timeDays=t;
T.config=cfg;

T.forcing=struct;
T.forcing.baselineDeltaDInterval_permil=d0;
T.forcing.plusDeltaDInterval_permil=dp;
T.forcing.minusDeltaDInterval_permil=dm;
T.forcing.pulseIntervalIndex=pulseIdx;
T.forcing.pulseStart_days=pulseStart;
T.forcing.pulseEnd_days=pulseEnd;
T.forcing.pulseDuration_days=pulseEnd-pulseStart;
T.forcing.pulseCentroid_days=inputCentroid;
T.forcing.inputOddAmplitude_permil=inputAmplitude;
T.forcing.inputOddArea_permil_days=inputArea;
T.forcing.inputEvenMax_permil=max(abs(inputEven));

T.response=response;
T.stateResponseSummary=stateResponseSummary;
T.causalTransferSummary=causalTransferSummary;

T.soil=struct;
T.soil.oddResponse_permil=oddSoil;
T.soil.evenResponse_permil=evenSoil;
T.soil.cellSummary=soilCellSummary;
T.soil.zoneResponse=soilZoneResponse;
T.soil.zoneSummary=soilZoneSummary;

T.wax=struct;
T.wax.response=waxResponse;
T.wax.homologSummary=waxHomologSummary;

T.closureSummary=closureSummary;

T.interpretation={ ...
    'Differential anomalies isolate isotope-forcing propagation because hydraulics and meteorology are identical among runs.', ...
    'Odd response is the mirror-symmetric tracer signal; even response measures nonlinear/asymmetric behavior.', ...
    'Peak transfer fraction is output peak isotope anomaly divided by the 100-per-mil precipitation pulse amplitude.', ...
    'State-response area gain compares time-integrated state anomaly with time-integrated precipitation forcing.', ...
    'Centroid lag is the timing shift of the positive state anomaly relative to the precipitation pulse centroid.', ...
    'Tail area fraction quantifies persistence after the precipitation pulse ends.', ...
    'Evaporation site is a diagnostic branch of the leaf system, not a causal through-flow stage between bulk leaf water and intracellular water.', ...
    'A substantial final residual flags truncation of the response by the available simulation window.'};

if cfg.makeFigures
    T.figureFiles=plot_full_model_precip_tracer_experiment( ...
        Mbase,Bbase,Mplus,Bplus,Mminus,Bminus,T,cfg);
else
    T.figureFiles={};
end

%% Console
fprintf('\n');
fprintf('============================================================\n');
fprintf('PAIRED PRECIPITATION-ISOTOPE TRACER ANALYSIS V1\n');
fprintf('============================================================\n');
fprintf('Pulse intervals: %d\n',numel(pulseIdx));
fprintf('Pulse start/end: %.6f -> %.6f days\n',pulseStart,pulseEnd);
fprintf('Pulse centroid:  %.6f days\n',inputCentroid);
fprintf('Positive pulse anomaly: +%.6f per mil\n',max(inputPlus));
fprintf('Negative pulse anomaly: %.6f per mil\n',min(inputMinus));
fprintf('Odd input amplitude: %.6f per mil\n',inputAmplitude);
fprintf('------------------------------------------------------------\n');
fprintf('STATE RESPONSE SUMMARY\n');
disp(stateResponseSummary);
fprintf('\nCAUSAL TRANSFER SUMMARY\n');
disp(causalTransferSummary);
fprintf('\nSOIL ZONE SUMMARY\n');
disp(soilZoneSummary);
fprintf('\nWAX HOMOLOG SUMMARY\n');
disp(waxHomologSummary);
fprintf('\nCLOSURE SUMMARY\n');
disp(closureSummary);
fprintf('============================================================\n');
end

%% =========================================================================
function cfg=setDefaultsLocal(cfg)
if ~isfield(cfg,'onsetFraction') || isempty(cfg.onsetFraction)
    cfg.onsetFraction=0.05;
end
if ~isfield(cfg,'truncationFinalFraction') || isempty(cfg.truncationFinalFraction)
    cfg.truncationFinalFraction=0.05;
end
if ~isfield(cfg,'makeFigures') || isempty(cfg.makeFigures)
    cfg.makeFigures=false;
end
if ~isfield(cfg,'outputDirectory') || isempty(cfg.outputDirectory)
    cfg.outputDirectory='full_model_tracer_output';
end
if ~isfield(cfg,'figureVisible') || isempty(cfg.figureVisible)
    cfg.figureVisible='off';
end
end

%% =========================================================================
function R=responseMetricsLocal( ...
    t,p,m,odd,even,inputAmplitude,inputArea, ...
    pulseStart,pulseEnd,inputCentroid,cfg)

t=t(:);
p=p(:);
m=m(:);
odd=odd(:);
even=even(:);

plusOriented=p;
minusOriented=-m;

positive=max(odd,0);
negative=max(-odd,0);

[peak,ix]=max(odd);
peakPlus=max(plusOriented);
peakMinus=max(minusOriented);

if ~isfinite(peak) || peak<=0
    peakTime=NaN;
    onsetTime=NaN;
else
    peakTime=t(ix);
    threshold=cfg.onsetFraction*peak;
    io=find(odd>=threshold,1,'first');
    if isempty(io)
        onsetTime=NaN;
    else
        onsetTime=t(io);
    end
end

positiveArea=trapz(t,positive);
negativeArea=trapz(t,negative);
netArea=trapz(t,odd);

if positiveArea>0
    centroid=trapz(t,t.*positive)./positiveArea;
else
    centroid=NaN;
end

tailMask=t>=pulseEnd;
if sum(tailMask)>=2 && positiveArea>0
    tailArea=trapz(t(tailMask),positive(tailMask));
    tailFraction=tailArea/positiveArea;
else
    tailFraction=NaN;
end

if peak>0
    finalResidual=odd(end)/peak;
else
    finalResidual=NaN;
end

if positiveArea>0
    negativeLobeFraction=negativeArea/positiveArea;
else
    negativeLobeFraction=NaN;
end

symMax=max(abs(even),[],'omitnan')/inputAmplitude;
symRMS=sqrt(mean(even.^2,'omitnan'))/inputAmplitude;
corrPM=correlationLocal(plusOriented,minusOriented);

R=struct;
R.peakPlusOriented_permil=peakPlus;
R.peakMinusOriented_permil=peakMinus;
R.peakOddResponse_permil=peak;
R.peakTransferFractionOfRain=peak/inputAmplitude;
R.peakTime_days=peakTime;
R.peakLagFromPulseCentroid_days=peakTime-inputCentroid;
R.onsetTime05_days=onsetTime;
R.onsetLagFromPulseStart_days=onsetTime-pulseStart;
R.centroidTime_days=centroid;
R.centroidLagFromPulse_days=centroid-inputCentroid;
R.positiveArea_permil_days=positiveArea;
R.netArea_permil_days=netArea;
R.areaGainVsRain=positiveArea/inputArea;
R.netAreaGainVsRain=netArea/inputArea;
R.tailAreaFractionAfterPulse=tailFraction;
R.finalResidualFractionOfPeak=finalResidual;
R.negativeLobeAreaFraction=negativeLobeFraction;
R.symmetryMaxResidualFractionOfInput=symMax;
R.symmetryRMSResidualFractionOfInput=symRMS;
R.plusVsOrientedMinusCorrelation=corrPM;
R.windowLikelyTruncated= ...
    isfinite(finalResidual) && abs(finalResidual)>cfg.truncationFinalFraction;
end

%% =========================================================================
function r=correlationLocal(x,y)
x=x(:); y=y(:);
ok=isfinite(x)&isfinite(y);
x=x(ok); y=y(ok);

if numel(x)<3 || std(x,0)<1e-12 || std(y,0)<1e-12
    r=NaN;
    return;
end

C=corrcoef(x,y);
r=C(1,2);
end

%% =========================================================================
function d=zoneBulkDeltaLocal(M,zoneNumber)

zone=M.grid.zoneIndex(:);
cells=zone==zoneNumber;

W=M.soilIsotopes.waterAmount_mol_m2(cells,:);
x=M.soilIsotopes.deuteriumAtomFraction(cells,:);

xBulk=sum(W.*x,1)./sum(W,1);
R=xBulk./(1-xBulk);

d=(1000.*(R./M.config.Rstandard_DH-1)).';
end

%% =========================================================================
function validateRunDimensionsLocal(M0,Mp,Mm,nT)

runs={M0,Mp,Mm};
for i=1:numel(runs)
    if numel(runs{i}.timeDays)~=nT
        error('Tracer run time dimensions differ.');
    end
end
end
