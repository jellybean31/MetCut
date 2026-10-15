function T = wax_timing_validation_targets(workbookPath)
%WAX_TIMING_VALIDATION_TARGETS
% Assemble timing/turnover targets for standing plant-wax validation.
%
% Workbook:
%   Tipple2013: ~2-week formation window; wax essentially fixed ~19 weeks
%   Sachse2010: ~20 per mil wax seasonal change vs ~40 per mil water variation;
%               variation occurs among leaf generations
%
% See:
%   Gamarra & Kahmen 2017: mature C3 grass secondary synthesis 0.09-1.09 % d^-1
%   Gao et al. 2012: C27-C31 recycling time 71-128 d in Phleum pratense

if nargin<1 || isempty(workbookPath)
    projectRoot=fileparts(fileparts(mfilename('fullpath')));
    workbookPath=fullfile(projectRoot,'data','reference', ...
        'leaf_wax_biosynthesis_validation_dataset_v2.xlsx');
end

if ~isfile(workbookPath)
    error('Validation workbook not found: %s',workbookPath);
end

L=readtable(workbookPath,'Sheet','Literature_Targets', ...
    'VariableNamingRule','preserve');

T=struct;
T.workbookPath=workbookPath;
T.workbookRows=L;

% Pull the workbook timing rows exactly by Study_ID.
tippleRow=L(strcmp(string(L.Study_ID),"Tipple2013"),:);
sachseRow=L(strcmp(string(L.Study_ID),"Sachse2010"),:);

if isempty(tippleRow)
    error('Tipple2013 timing target not found in Literature_Targets.');
end

if isempty(sachseRow)
    error('Sachse2010 timing target not found in Literature_Targets.');
end

T.workbook.Tipple2013=tippleRow;
T.workbook.Sachse2010=sachseRow;

% Numeric operationalization of the workbook timing statements.
T.operational.Tipple2013.formationWindowDays=14;
T.operational.Tipple2013.postFormationStableDays=19*7;
T.operational.Tipple2013.description= ...
    ['Deciduous broadleaf finite production window: approximately 2 weeks, ' ...
     'followed by approximately 19 weeks of essentially fixed plant-wax deltaD.'];

T.operational.Sachse2010.waxSeasonalAmplitude_permil=20;
T.operational.Sachse2010.waterSeasonalAmplitude_permil=40;
T.operational.Sachse2010.description= ...
    ['Barley plant wax varies among leaf generations produced at different ' ...
     'times, with little within-generation change after formation.'];

% External literature targets, kept explicitly separate from workbook rows.
T.external.Gamarra2017.study='Gamarra & Kahmen 2017';
T.external.Gamarra2017.DOI='10.1002/rcm.7770';
T.external.Gamarra2017.matureSecondarySynthesisFractionPerDay=[0.0009 0.0109];
T.external.Gamarra2017.description= ...
    'Six mature C3 grass species; secondary n-alkane synthesis 0.09-1.09% d^-1.';

T.external.Gao2012.study='Gao et al. 2012';
T.external.Gao2012.DOI='10.1002/rcm.5313';
T.external.Gao2012.longChainRecycleTimeDays=[71 128];
T.external.Gao2012.approxFirstOrderRatePerDay=1./[128 71];
T.external.Gao2012.description= ...
    ['Phleum pratense D-labeling; C27-C31 recycling time 71-128 d. ' ...
     '1/tau conversion is used only as an approximate first-order mapping.'];

end
