function P = inference_prior_registry()
%INFERENCE_PRIOR_REGISTRY
%
% Uncertain model dimensions and prior status.
%
% PRIOR CLASSES
% -------------
% development
%   Broad development/sensitivity-analysis prior or envelope.
%
% development_test
%   Working distribution used for synthetic/inference tests only.
%
% empirical
%   Numerical prior directly supported by calibration data.
%
% conditional
%   Numerical prior conditioned on external PFT/taxon/trait/environment data.
%
% absent
%   No numerical prior currently assigned.
%
% INFERENCE ROLES
% ---------------
% infer
% marginalize
% infer_or_marginalize
% marginalize_or_condition
% condition
% optional_zero_default
% optional_off_default
% optional_dynamic


rows = { ...
% id                     parameterIDs                  priorClass         numericalStatus       inferenceRole                 priority     mandatory  optional  conditioningVariables                                  evidenceSource
'BIO_MODE12',            {'B01','B02','B03','B04','B05','B06','B07'}, 'development',      'available',           'infer_or_marginalize',          'HIGH',     true,      false,    {'observation_set'},                                     {'K3.2E','M1','M2'};
'B08_W01',               {'B08','W01'},                              'development',      'available',           'marginalize',                   'LOW',      true,      false,    {'taxon','leaf_development_state'},                      {'J1','M3','M4b'};
'I01_I02',               {'I01','I02'},                              'development',      'available',           'marginalize_or_condition',      'MEDIUM',   true,      false,    {'plant_size','PFT','growth_form'},                      {'K3.2C','M3','M4b'};
'I07_I08',               {'I07','I08'},                              'development',      'available',           'marginalize',                   'LOW',      true,      false,    {'PFT','leaf_cell_anatomy'},                             {'J1','M3','M4b'};
'ROOTING_ARCH',          {'ROOTING_ARCH'},                           'conditional',      'available',           'condition',                     'HIGH',     true,      false,    {'PET','MAP_over_PET','growth_form'},                    {'L2b.1','L2f','M4b','M6d','M6g','M6h','M6i','M6j','M6k'};
'LEAF_ANATOMY',          {'LEAF_ANATOMY'},                           'absent',           'none',                'condition',                     'HIGH',     true,      false,    {'PFT','taxon','leaf_thickness','mesophyll_structure','vein_density'},       {'K3.2A','M4b'};
'W02',                   {'W02'},                                    'development',      'available',           'marginalize_or_condition',      'MEDIUM',   true,      false,    {'taxon','leaf_age','maturity','turnover_state'},        {'J1','K3','M4a'};
'ROOT_FRAC',             {'ROOT_FRAC'},                              'absent',           'none',                'optional_zero_default',         'OPTIONAL', false,     true,     {'species','root_transport_evidence'},                   {'K3.2D','M4b'};
'SOIL_EVAP_ENRICH',      {'SOIL_EVAP_ENRICH'},                       'absent',           'none',                'optional_off_default',          'OPTIONAL', false,     true,     {'aridity','soil_evaporation','soil_depth','event_history'},               {'K3.2D','L2','M4b'};
'DYNAMIC_RESERVE_C',     {'DYNAMIC_RESERVE_C'},                      'absent',           'none',                'optional_dynamic',              'OPTIONAL', false,     true,     {'leaf_development_state','pulse_label','metabolic_state'},                {'K3.2F','M4b'};
};

n=size(rows,1);

P=repmat(struct( ...
'id',"", ...
'parameterIDs',strings(0,1), ...
'priorClass',"", ...
'numericalPriorStatus',"", ...
'inferenceRole',"", ...
'priority',"", ...
'mandatoryInArchitecture',false, ...
'optionalProcess',false, ...
'conditioningVariables',strings(0,1), ...
'marginalizationPermitted',false, ...
'freeInferencePermitted',false, ...
'requiresExternalConditioning',false, ...
'defaultEnabled',true, ...
'evidenceSource',strings(0,1), ...
'notes',""),n,1);

for i=1:n
P(i).id=string(rows{i,1});
P(i).parameterIDs=string(rows{i,2}(:));
P(i).priorClass=string(rows{i,3});
P(i).numericalPriorStatus=string(rows{i,4});
P(i).inferenceRole=string(rows{i,5});
P(i).priority=string(rows{i,6});
P(i).mandatoryInArchitecture=logical(rows{i,7});
P(i).optionalProcess=logical(rows{i,8});
P(i).conditioningVariables=string(rows{i,9}(:));
P(i).evidenceSource=string(rows{i,10}(:));

role=P(i).inferenceRole;

P(i).marginalizationPermitted = any(role==[ ...
    "marginalize","infer_or_marginalize","marginalize_or_condition"]);

P(i).freeInferencePermitted = any(role==[ ...
    "infer","infer_or_marginalize"]);

P(i).requiresExternalConditioning = any(role==[ ...
    "condition","marginalize_or_condition"]);

if P(i).optionalProcess
    P(i).defaultEnabled=false;
else
    P(i).defaultEnabled=true;
end
end

% Row-specific interpretation notes.
setNote('BIO_MODE12', ...
'Mode 1 is the routine C29 biochemical coordinate; Mode 2 is added for multi-homologue inference. Full B01-B07 remain mechanistic forward parameters.');
setNote('B08_W01', ...
'M3 shows negligible degradation when marginalized under current development priors.');
setNote('I01_I02', ...
'Marginalize by default for mature standing wax; condition on plant size/PFT where independent information exists.');
setNote('I07_I08', ...
'Marginalize under current development priors for routine mature standing-wax inference.');
setNote('ROOTING_ARCH', ...
'Validated conditional logistic-normal simplex mapping from RPGE root profiles using PET, MAP/PET, and growth form. Robust Huber ALR mean with tail-calibrated covariance; sampling depth and below-0.9-m mass remain diagnostic metadata, not deployable predictors. Extrapolation is disabled outside the validated mapping domain.');
setNote('LEAF_ANATOMY', ...
'High-priority conditional dimension. No universal numerical PFT prior is currently assigned.');
setNote('W02', ...
'Taxon/leaf-age information may justify conditional narrowing; otherwise retain development prior.');
setNote('ROOT_FRAC', ...
'Zero default. Activate only with species/process evidence supporting nonzero apparent root uptake fractionation.');
setNote('SOIL_EVAP_ENRICH', ...
'Off by default. Activate only in applications where soil evaporative isotope enrichment is physically relevant.');
setNote('DYNAMIC_RESERVE_C', ...
'Off for mature endpoint inference; enable for developing leaves, pulse-label studies, or explicit dynamic metabolic questions.');

function setNote(id,note)
    j=find(string({P.id})==string(id),1);
    if isempty(j)
        error('Internal registry row %s not found.',id);
    end
    P(j).notes=string(note);
end
end
