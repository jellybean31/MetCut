function R = test_METCUT_manuscript_wax_timing_targets()
%TEST_METCUT_MANUSCRIPT_WAX_TIMING_TARGETS
%
% Regression for standing plant-wax timing/turnover.


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

workbookPath=fullfile(refDir, ...
    'leaf_wax_biosynthesis_validation_dataset_v2.xlsx');

assert(isfile(workbookPath), ...
    ['Package reference workbook missing. Place ' ...
     'leaf_wax_biosynthesis_validation_dataset_v2.xlsx in data/reference/.']);

T=wax_timing_validation_targets();

assert(isfile(T.workbookPath), ...
    'wax_timing_validation_targets returned a nonexistent workbook path.');

assert(strcmp(canonicalLocal(T.workbookPath),canonicalLocal(workbookPath)), ...
    'Timing-target helper did not resolve the package-local reference workbook.');

% -------------------------------------------------------------------------
% Workbook-derived rows.
% -------------------------------------------------------------------------
assert(height(T.workbook.Tipple2013)==1, ...
    'Expected exactly one Tipple2013 timing row.');
assert(height(T.workbook.Sachse2010)==1, ...
    'Expected exactly one Sachse2010 timing row.');

assert(string(T.workbook.Tipple2013.Study_ID)=="Tipple2013", ...
    'Tipple2013 Study_ID regression.');
assert(string(T.workbook.Sachse2010.Study_ID)=="Sachse2010", ...
    'Sachse2010 Study_ID regression.');

assert(ismember('Type',T.workbook.Tipple2013.Properties.VariableNames), ...
    'Literature_Targets release is missing the Type column.');
assert(string(T.workbook.Tipple2013.Type)=="Timing", ...
    'Tipple2013 timing type regression.');
assert(ismember('Type',T.workbook.Sachse2010.Properties.VariableNames), ...
    'Literature_Targets release is missing the Type column.');
assert(string(T.workbook.Sachse2010.Type)=="Seasonality", ...
    'Sachse2010 seasonality type regression.');

tippleObservation=string(T.workbook.Tipple2013.Observation);
tippleValue=string(T.workbook.Tipple2013.Value);
tippleRange=string(T.workbook.Tipple2013.Uncertainty_or_Range);

sachseObservation=string(T.workbook.Sachse2010.Observation);
sachseValue=string(T.workbook.Sachse2010.Value);
sachseRange=string(T.workbook.Sachse2010.Uncertainty_or_Range);

assert(contains(lower(tippleObservation),"leaf flush"), ...
    'Tipple2013 workbook observation regression.');
assert(contains(tippleValue,"2-week"), ...
    'Tipple2013 workbook formation-window text regression.');
assert(contains(tippleRange,"19 weeks"), ...
    'Tipple2013 workbook stability-window text regression.');

assert(contains(sachseValue,"20"), ...
    'Sachse2010 workbook wax-amplitude text regression.');
assert(contains(sachseRange,"40"), ...
    'Sachse2010 workbook water-amplitude text regression.');
assert(contains(lower(sachseObservation),"season"), ...
    'Sachse2010 workbook observation regression.');

% -------------------------------------------------------------------------
% Operationalized manuscript targets.
% -------------------------------------------------------------------------
assert(T.operational.Tipple2013.formationWindowDays==14, ...
    'Tipple2013 14-day formation-window regression.');
assert(T.operational.Tipple2013.postFormationStableDays==133, ...
    'Tipple2013 19-week stability-window regression.');

assert(T.operational.Sachse2010.waxSeasonalAmplitude_permil==20, ...
    'Sachse2010 plant-wax amplitude regression.');
assert(T.operational.Sachse2010.waterSeasonalAmplitude_permil==40, ...
    'Sachse2010 water-amplitude regression.');

% -------------------------------------------------------------------------
% External timing/turnover constraints remain explicitly separate.
% -------------------------------------------------------------------------
gamarra=double(T.external.Gamarra2017.matureSecondarySynthesisFractionPerDay);
gaoTau=double(T.external.Gao2012.longChainRecycleTimeDays);
gaoRate=double(T.external.Gao2012.approxFirstOrderRatePerDay);

assert(max(abs(gamarra-[0.0009 0.0109]))<=1e-15, ...
    'Gamarra & Kahmen secondary-synthesis range regression.');
assert(isequal(gaoTau,[71 128]), ...
    'Gao long-chain recycling-time regression.');
assert(max(abs(gaoRate-(1./[128 71])))<=1e-15, ...
    'Gao approximate first-order-rate mapping regression.');

assert(strcmp(T.external.Gamarra2017.DOI,'10.1002/rcm.7770'), ...
    'Gamarra DOI regression.');
assert(strcmp(T.external.Gao2012.DOI,'10.1002/rcm.5313'), ...
    'Gao DOI regression.');

fprintf('\n');
fprintf('============================================================\n');
fprintf('METCUT STANDING-WAX TIMING-TARGET REGRESSION\n');
fprintf('============================================================\n');
fprintf('Workbook:                         %s\n',T.workbookPath);
fprintf('Tipple formation window:          %.6g d\n', ...
    T.operational.Tipple2013.formationWindowDays);
fprintf('Tipple post-formation stability:  %.6g d\n', ...
    T.operational.Tipple2013.postFormationStableDays);
fprintf('Sachse wax/water amplitude:       %.6g / %.6g permil\n', ...
    T.operational.Sachse2010.waxSeasonalAmplitude_permil, ...
    T.operational.Sachse2010.waterSeasonalAmplitude_permil);
fprintf('Gamarra secondary synthesis:      %s d^-1\n',mat2str(gamarra,12));
fprintf('Gao recycling time:               %s d\n',mat2str(gaoTau));
fprintf('PASS\n');
fprintf('============================================================\n');

R=struct;
R.pass=true;
R.workbookPath=T.workbookPath;
R.TippleFormationWindowDays=T.operational.Tipple2013.formationWindowDays;
R.TippleStableDays=T.operational.Tipple2013.postFormationStableDays;
R.SachseWaxAmplitude_permil=T.operational.Sachse2010.waxSeasonalAmplitude_permil;
R.SachseWaterAmplitude_permil=T.operational.Sachse2010.waterSeasonalAmplitude_permil;
R.GamarraSecondarySynthesisFractionPerDay=gamarra;
R.GaoRecycleTimeDays=gaoTau;
R.GaoApproxRatePerDay=gaoRate;
end

% =========================================================================
function p=canonicalLocal(p)
p=char(java.io.File(p).getCanonicalPath());
end
