# METCUT v1 Tutorial

## Purpose

METCUT is a coupled process model for propagating hydrogen-isotope signals from water inputs through soil, roots, xylem, leaf water, intracellular water pools, plant-wax biochemistry, and standing plant-wax accumulation.

This tutorial is designed for a user who has downloaded the release package and wants to:

1. verify that MATLAB can see the release correctly;
2. run individual model modules;
3. inspect the intermediate states produced by the model;
4. run the complete forward model;
5. perturb selected inputs;
6. run standing-wax and inverse examples;
7. reproduce validation tests.

---

# 1. Expected layout

The examples below assume a release directory with this general structure:

METCUT_v1/
├── data/
│   └── reference/
├── functions/
├── tests/
│   ├── unit/
│   ├── integration/
│   └── manuscript/
└── METCUT_v1_Tutorial.md


The runtime uses the payloads under 'data/reference/'.

The large development-stage 'data/input/' tree is not required by the current validated runtime.

---

# 2. MATLAB setup

Open MATLAB and change the current folder to the top-level 'METCUT_v1' directory.

For example:

```
cd('/path/to/METCUT_v1')
```

Then add the production functions:

```
clear functions
rehash

root = pwd;
fnDir = fullfile(root,'functions');

addpath(fnDir,'-begin');
```

Check several important release interfaces:

```
which run_full_model_from_hydraulic_checkpoint
which effective7_node_wax
which intracellular_water_dynamic_exactD
which wax_accumulation_dynamic_exactD
which build_inverse_engine_spec
which inverse_runtime_dispatcher
```

All should resolve under:

```
METCUT_v1/functions/
```

---

# 3. Basic model architecture 

The main process sequence is:

```
precipitation / initial soil water
        ↓
soil water balance and isotope transport
        ↓
root-source weighting
        ↓
plant hydraulics
        ↓
xylem isotope transport and storage
        ↓
leaf evaporative enrichment
        ↓
intracellular cytosol / plastid water
        ↓
biochemical hydrogen-source partitioning
        ↓
instantaneous C27 / C29 / C31 / C33 plant wax
        ↓
standing-wax production, replacement, loss, and memory
```


---

# 4. Recommended way to learn the modules

Some METCUT components are simple numerical primitives and can be called directly.

Other components require large hydraulic, meteorological, anatomical, or biochemical state structures. For those modules, the runnable examples are the tests in 'tests/unit/' and 'tests/integration/'. Those tests construct inputs and check exact-D closure and numerical behavior.

---

# 5. Exact-D numerical primitives

These are the lowest-level isotope utilities used throughout METCUT.

Important functions include:

```
delta2ratio
ratio2delta
exact_epsilon_permil
effective_state_exactD_mix
```

Example:

```
addpath(fullfile(root,'tests','unit'),'-begin');

R = test_METCUT_effective_state_exactD_mix;
```

This checks exact source and water endpoints and monotonic interpolation.

You can also inspect the function documentation:

```
help effective_state_exactD_mix
help delta2ratio
help ratio2delta
```

Remove the unit-test directory afterward if desired:

```
rmpath(fullfile(root,'tests','unit'));
```

---

# 6. Soil water and soil-isotope module

The retained soil architecture includes:

```
build_richards_soil_grid
soil_matric_flux_potential
soil_relative_conductivity
soil_root_layer_flux
soil_water_balance_richards_1d_prescribed_root
soil_isotope_transport_richards
```

The soil solver is coupled to plant hydraulics in the validated release, so the runner is the hydraulic integration regression.

```
addpath(fullfile(root,'tests','integration'),'-begin');

Rhyd = test_METCUT_hydraulic_solver_regression;
```

This runs the hydraulic solver. The hydraulic reference used by the quicker examples below is:

```
hydPath = fullfile(root,'data','reference', ...
    'hydraulic_checkpoint_76cell_8day_hourly_corrected_mass_conservative_storagefixed.mat');
```

---

# 7. Plant-hydraulic module

Important retained hydraulic functions include:

```
solve_plant_hydraulics
couple_richards_soil_plant_hydraulics
hydraulic_vulnerability
xylem_storage_pressure_curve
xylem_storage_hysteretic
transpirationDemand
```

Because the coupled solver requires a complete soil/root/meteorological configuration, use the regression as the first runnable hydraulic example:

```
addpath(fullfile(root,'tests','integration'),'-begin');

Rhyd = test_METCUT_hydraulic_solver_regression;
```

For source inspection:

```
help solve_plant_hydraulics
help couple_richards_soil_plant_hydraulics
```


---

# 8. Xylem isotope transport

The isotope-transport primitive is:

```
xylem_transport_dynamic_exactD
```

This represents isotope transport through the mobile xylem state, including the storage/recharge/release history represented by the model.

The simplest way to inspect the resulting xylem trajectory is to run the complete reference model:

```
[M,B] = run_full_model_from_hydraulic_checkpoint(hydPath);

xylemD = M.xylemIsotopes.leafBoundary.deltaDWater;

plot(M.timeDays,xylemD)
xlabel('Time (days)')
ylabel('Xylem-delivered \deltaD (‰)')
title('METCUT xylem isotope trajectory')
```

The complete model also exposes the exact-D xylem closure diagnostic:

```
maxXylemClosure = max(abs( ...
    M.xylemIsotopes.deuteriumInventoryClosureResidual_molD_m2(:)), ...
    [],'omitnan');

fprintf('Maximum xylem D closure residual = %.3e mol D m^-2\n', ...
    maxXylemClosure);
```

---

# 9. Leaf energy balance and atmospheric exchange

Retained leaf-surface and energy-balance functions include:

```
mixed_surface_heat
mixed_surface_scalar
leaf_energy_balance_residual
leaf_energy_balance_root_solve
surface_to_scalar_gb_forward
air_transport_diffusivities
water_vapor_diffusivity_air_massman1998
```

These functions support the leaf-water module but require coupled meteorological and surface-state inputs.

For first use, inspect them through the full model and integrations:

```
[M,B] = run_full_model_from_hydraulic_checkpoint(hydPath);

fieldnames(M.config)
```

The hydraulic and full-model integration tests exercise the coupled exchange configuration.

---

# 10. Leaf evaporative-enrichment module

The main leaf-water function is:

```
leaf_hydraulics_dynamic_spatial_symplast
```

In the coupled architecture, this module receives:

```
time
xylem-delivered isotope state
hydraulic state
meteorological / boundary-layer state
leaf isotope traits
```

Its output is stored in the full model as:

```
M.leafIsotopes
```

Run the reference model and inspect the bulk symplast and evaporation-site signals:

```
[M,B] = run_full_model_from_hydraulic_checkpoint(hydPath);

t = M.timeDays(:);

bulkLeaf = M.leafIsotopes.deltaDSymplastBulk(:);

plot(t,bulkLeaf)
xlabel('Time (days)')
ylabel('Bulk symplast \deltaD (‰)')
title('Leaf-water isotope response')
```

The exact-D leaf closure can be inspected with:

```
maxLeafClosure = max(abs( ...
    M.leafIsotopes.DInventoryClosureResidual_molD_m2(:)), ...
    [],'omitnan');

fprintf('Maximum leaf D closure residual = %.3e mol D m^-2\n', ...
    maxLeafClosure);
```

---

# 11. Intracellular-water module

The intracellular module is:

```
intracellular_water_dynamic_exactD
```

It propagates water and D inventories through the retained intracellular water pools, including cytosol and plastid states.

A numerical validation is:

```
addpath(fullfile(root,'tests','unit'),'-begin');

Rcell = test_METCUT_intracellular_water_exactD;
```

To inspect the intracellular trajectories produced by the complete model:

```
[M,B] = run_full_model_from_hydraulic_checkpoint(hydPath);

t = M.timeDays(:);
cy = M.intracellular.deltaD.cytosol(:);
pl = M.intracellular.deltaD.plastid(:);

plot(t,cy,'DisplayName','Cytosol')
hold on
plot(t,pl,'DisplayName','Plastid')
hold off

xlabel('Time (days)')
ylabel('\deltaD (‰)')
legend('Location','best')
title('Intracellular water pools')
```

The exact-D system closure is:

```
maxCellClosure = max(abs( ...
    M.intracellular.wholeSystemDClosureResidual_molD_m2(:)), ...
    [],'omitnan');

fprintf('Maximum intracellular D closure residual = %.3e mol D m^-2\n', ...
    maxCellClosure);
```

---

# 12. Biochemical hydrogen-incorporation module

The biochemical architecture includes:

```
alkane_biosynthesis_h_stoichiometry_reaction_explicit
wax_hydrogen_source_pools_exactD
wax_combine_hydrogen_classes_exactD
wax_apply_reaction_incorporation_exactD
effective7_node_wax
effective7_structural_basis
biochemical_latent_mode12_transform
```

The core exact-D reaction-incorporation step has a unit example:

```
addpath(fullfile(root,'tests','unit'),'-begin');

Rbio = test_METCUT_reaction_incorporation_exactD;
```

The full model stores the instantaneous plant-wax isotope trajectories under:

```
M.biochemistry.wax
```

Run the complete reference model:

```
[M,B] = run_full_model_from_hydraulic_checkpoint(hydPath);

W = double(M.biochemistry.wax.deltaDMatrix);
```

The baseline contains four homologues. If explicit carbon numbers are present, use them. Otherwise the validated fallback row order is C27, C29, C31, C33:

```
if isfield(M.biochemistry.wax,'alkaneCarbonNumber')
    C = double(M.biochemistry.wax.alkaneCarbonNumber(:).');
else
    C = [27 29 31 33];
end

if size(W,1) ~= numel(C) && size(W,2) == numel(C)
    W = W.';
end
```

Plot all four instantaneous plant-wax trajectories:

```
plot(M.timeDays,W.')
xlabel('Time (days)')
ylabel('Instantaneous plant-wax \deltaD (‰)')
legend(compose('C%d',C),'Location','best')
title('Instantaneous plant-wax isotope trajectories')
```

---

# 13. Standing plant-wax accumulation

The standing-wax functions include:

```
wax_accumulation_default_traits
wax_accumulation_production_schedule
wax_accumulation_dynamic_exactD
attach_wax_accumulation_to_model
```

A exact-D validation is:

```
addpath(fullfile(root,'tests','unit'),'-begin');

Rwax = test_METCUT_standing_wax_accumulation_exactD;
```

The coupled integration example is:

```
addpath(fullfile(root,'tests','integration'),'-begin');

RwaxInt = test_METCUT_standing_wax_integration;
```

A typical accumulation-trait configuration begins with:

```
accumTraits = wax_accumulation_default_traits();

accumTraits.productionMode = "developmentPlusMaintenance";
accumTraits.developmentWindowDays = 14;
accumTraits.matureReplacementFractionPerDay = sqrt(0.0009*0.0109);
accumTraits.lossRatePerDay = sqrt(0.0009*0.0109);
```

These values are useful as a worked example.

For experiments involving standing wax, the integration test is the best template because it uses the correct model stage and production schedule consistently.

---

# 14. Full forward model: quickest complete run

This is the main starting point for most users.

```
clear functions
rehash

root = pwd;
addpath(fullfile(root,'functions'),'-begin');

hydPath = fullfile(root,'data','reference', ...
    'hydraulic_checkpoint_76cell_8day_hourly_corrected_mass_conservative_storagefixed.mat');

assert(isfile(hydPath),'Hydraulic reference is missing.');

[M,B] = run_full_model_from_hydraulic_checkpoint(hydPath);
```

The returned model structure contains the major stages:

```
fieldnames(M)
```

Important outputs include:

```
M.soilIsotopes
M.rootBoundary
M.xylemIsotopes
M.leafIsotopes
M.intracellular
M.biochemistry
```

A useful quick summary is:

```matlab
fprintf('Final root-source dD:     %.6f per mil\n', ...
    M.rootBoundary.deltaDNode_permil(end));

fprintf('Final xylem dD:           %.6f per mil\n', ...
    M.xylemIsotopes.leafBoundary.deltaDWater(end));

fprintf('Final bulk-leaf dD:       %.6f per mil\n', ...
    M.leafIsotopes.deltaDSymplastBulk(end));

fprintf('Final cytosol dD:         %.6f per mil\n', ...
    M.intracellular.deltaD.cytosol(end));

fprintf('Final plastid dD:         %.6f per mil\n', ...
    M.intracellular.deltaD.plastid(end));
```

Extract final C27/C29/C31/C33 values:

```
W = double(M.biochemistry.wax.deltaDMatrix);

if isfield(M.biochemistry.wax,'alkaneCarbonNumber')
    C = double(M.biochemistry.wax.alkaneCarbonNumber(:).');
else
    C = [27 29 31 33];
end

if size(W,1) ~= numel(C) && size(W,2) == numel(C)
    W = W.';
end

finalWax = NaN(size(C));

for i = 1:numel(C)
    row = W(i,:);
    j = find(isfinite(row),1,'last');
    finalWax(i) = row(j);
end

disp(table(C(:),finalWax(:), ...
    'VariableNames',{'AlkaneC','FinalDeltaD_permil'}))
```

The regression is:

```
C27  -146.557420 ‰
C29  -146.344758 ‰
C31  -146.158680 ‰
C33  -145.994493 ‰
```

---

# 15. Parameterized full-model run

The release contains a parameterized sibling:

```text
run_full_model_from_hydraulic_checkpoint_parameterized
```

The integration test provides an example parameter structure:

```
P = struct;

P.initialSoilDeltaDByZone_permil = [-60; -75; -90];
P.precipDeltaD_permil = -50;

P.leafAnatomyConditioning = struct( ...
    'interveinalDistance_m',0.0005, ...
    'mesophyllThickness_m',0.00015, ...
    'effectiveFlowAreaFraction',0.02, ...
    'flowFraction',1, ...
    'tortuosity',1.5, ...
    'lateralPathFraction',0.5, ...
    'verticalPathFraction',0.5);

P.biochemistryPolicy = struct;
P.biochemistryPolicy.parameterMode = 'coreMidpointDiagnostic';
P.biochemistryPolicy.routeFractions = struct('C16',0,'C18',1);
P.biochemistryPolicy.alkaneCarbonNumber = [27 29 31 33];

[Mp,Bp] = run_full_model_from_hydraulic_checkpoint_parameterized( ...
    hydPath,P);
```

The values above reproduce the reference defaults.

To perform a simple sensitivity experiment, modify one field at a time and compare the results.

For example, change precipitation δD:

```
P2 = P;
P2.precipDeltaD_permil = -70;

[M2,B2] = run_full_model_from_hydraulic_checkpoint_parameterized( ...
    hydPath,P2);
```

Compare the C29 response:

```
W1 = double(Mp.biochemistry.wax.deltaDMatrix);
W2 = double(M2.biochemistry.wax.deltaDMatrix);

if size(W1,1) ~= 4 && size(W1,2) == 4
    W1 = W1.';
end

if size(W2,1) ~= 4 && size(W2,2) == 4
    W2 = W2.';
end

fprintf('Reference final C29: %.6f per mil\n',W1(2,end));
fprintf('Changed final C29:   %.6f per mil\n',W2(2,end));
fprintf('Difference:          %.6f per mil\n',W2(2,end)-W1(2,end));
```

---

# 16. Inspecting process-by-process propagation

After a full run, a useful exercise is to plot one representative state from each stage.

```
[M,B] = run_full_model_from_hydraulic_checkpoint(hydPath);

t = M.timeDays(:);

rootD = M.rootBoundary.deltaDNode_permil(:);
xylemD = M.xylemIsotopes.leafBoundary.deltaDWater(:);
leafD = M.leafIsotopes.deltaDSymplastBulk(:);
cyD = M.intracellular.deltaD.cytosol(:);
plD = M.intracellular.deltaD.plastid(:);

W = double(M.biochemistry.wax.deltaDMatrix);

if size(W,1) ~= 4 && size(W,2) == 4
    W = W.';
end

c29 = W(2,:).';

plot(t,rootD,'DisplayName','Root source')
hold on
plot(t,xylemD,'DisplayName','Xylem delivered')
plot(t,leafD,'DisplayName','Leaf symplast')
plot(t,cyD,'DisplayName','Cytosol')
plot(t,plD,'DisplayName','Plastid')
plot(t,c29,'DisplayName','C29 instantaneous wax')
hold off

xlabel('Time (days)')
ylabel('\deltaD (‰)')
legend('Location','best')
title('METCUT process-by-process isotope propagation')
```


---

# 17. Precipitation-isotope tracer example

This test provides a complete perturbation example that propagates paired positive and negative precipitation-isotope pulses through the coupled model.

Run:

```
addpath(fullfile(root,'tests','manuscript'),'-begin');

Rtracer = test_METCUT_manuscript_precip_tracer;
```

This demonstrates propagation through:

```
soil
root source
xylem
leaf water
intracellular water
instantaneous C27/C29/C31/C33 plant wax
```

It is a good starting template for designing a time-varying isotope forcing experiment.

---

# 18. Source and temporal-memory example

To demonstrate source overrides and transport / standing-wax memory:

```
addpath(fullfile(root,'tests','manuscript'),'-begin');

Rmemory = test_METCUT_manuscript_source_memory;
```

This exercises workflows including:

```
rerun_full_model_isotope_source_overrides_with_wax_timing
rerun_full_model_temporal_memory_source_case
```

Use this test as the reference template before modifying residence-time or source-state assumptions.

---

# 19. Empirical timing targets

The timing helper reads the local workbook:

```
T = wax_timing_validation_targets();
```

The workbook is:

```
data/reference/leaf_wax_biosynthesis_validation_dataset_v2.xlsx
```

Inspect the returned structure:

```
fieldnames(T)
```

The regression is:

```
addpath(fullfile(root,'tests','manuscript'),'-begin');

Rtiming = test_METCUT_manuscript_wax_timing_targets;
```

The worked targets include:

```
Tipple formation window:          14 d
Tipple post-formation stability:  133 d
Sachse wax/water amplitude:       20 / 40 ‰
Gamarra secondary synthesis:      [0.0009 0.0109] d^-1
Gao recycling time:               [71 128] d
```


---

# 20. Inverse-model architecture

The inverse entry interfaces are:

```
build_inverse_engine_spec
inverse_runtime_dispatcher
```

Supporting components include:

```
inverse_configuration_template
resolve_inverse_configuration
apply_nuisance_transform
nuisance_transform_spec
gaussian_observation_likelihood
effective7_structural4_linear_observation_objective
```

The inverse engine expects a reference model, inverse configuration, and proposal state. Because those structures are conditional on the selected observation and conditioning design, a good first executable example is the full-chain recovery test:

```
addpath(fullfile(root,'tests','manuscript'),'-begin');

Rinv = test_METCUT_manuscript_inverse_recovery;
```

This performs a synthetic full-chain wax recovery experiment over:

```
I02
L1
L2
```

for both final and time-resolved plant-wax observations.


---

# 21. Running all validation examples

The following block runs the six unit tests, three integration tests, and four manuscript tests individually.


```
clear functions
rehash

root = pwd;

addpath(fullfile(root,'functions'),'-begin');
addpath(fullfile(root,'tests','unit'),'-begin');
addpath(fullfile(root,'tests','integration'),'-begin');
addpath(fullfile(root,'tests','manuscript'),'-begin');

%% Unit tests
test_METCUT_effective_state_exactD_mix;
test_METCUT_intracellular_water_exactD;
test_METCUT_reaction_incorporation_exactD;
test_METCUT_standing_wax_accumulation_exactD;
test_METCUT_subspace_principal_angles;
test_METCUT_synthetic_recovery_definitions;

%% Integration tests
test_METCUT_full_model_reference_regression;
test_METCUT_hydraulic_solver_regression;
test_METCUT_standing_wax_integration;

%% Manuscript / scientific regression tests
test_METCUT_manuscript_inverse_recovery;
test_METCUT_manuscript_precip_tracer;
test_METCUT_manuscript_source_memory;
test_METCUT_manuscript_wax_timing_targets;
```

Expected result:

```
Unit:        6 / 6 PASS
Integration: 3 / 3 PASS
Manuscript:  4 / 4 PASS
```

The hydraulic solver and inverse-recovery tests are the slowest of these examples.

---

# 22. A compact workflow

For a first session with METCUT, the recommended sequence is:

## Step 1 — run the complete model

```
root = pwd;
addpath(fullfile(root,'functions'),'-begin');

hydPath = fullfile(root,'data','reference', ...
    'hydraulic_checkpoint_76cell_8day_hourly_corrected_mass_conservative_storagefixed.mat');

[M,B] = run_full_model_from_hydraulic_checkpoint(hydPath);
```

## Step 2 — inspect the major model states

```
fieldnames(M)

M.soilIsotopes
M.rootBoundary
M.xylemIsotopes
M.leafIsotopes
M.intracellular
M.biochemistry
```

## Step 3 — plot process propagation

Use the example in Section 16.

## Step 4 — modify one parameterized input

Use the parameterized example in Section 15.

## Step 5 — run the precipitation tracer

```
addpath(fullfile(root,'tests','manuscript'),'-begin');
Rtracer = test_METCUT_manuscript_precip_tracer;
```

## Step 6 — run the standing-wax integration

```
addpath(fullfile(root,'tests','integration'),'-begin');
RwaxInt = test_METCUT_standing_wax_integration;
```

## Step 7 — run the inverse demonstration

```
addpath(fullfile(root,'tests','manuscript'),'-begin');
Rinv = test_METCUT_manuscript_inverse_recovery;
```

At that point you have exercised the main forward, temporal-memory, and inverse components of the model.


# 23. Useful model fields

After:

```
[M,B] = run_full_model_from_hydraulic_checkpoint(hydPath);
```

the following fields are useful for exploration:

```
M.timeDays

M.soilIsotopes
M.rootBoundary

M.hydraulics
M.xylemIsotopes

M.leafIsotopes
M.intracellular

M.biochemistry
M.biochemistry.wax

M.config
```

MATLAB's 'fieldnames' command is useful for exploring nested structures:

```
fieldnames(M)
fieldnames(M.xylemIsotopes)
fieldnames(M.leafIsotopes)
fieldnames(M.intracellular)
fieldnames(M.biochemistry)
fieldnames(M.biochemistry.wax)
```

---


# 24. Suggested tutorial exercises

## Exercise A — precipitation isotope sensitivity

1. Run the parameterized reference.
2. Change 'P.precipDeltaD_permil'.
3. Compare root, xylem, leaf, plastid, and C29 responses.
4. Determine which stages attenuate or delay the perturbation.

## Exercise B — initial soil isotope gradient

Change:

```
P.initialSoilDeltaDByZone_permil
```

while keeping precipitation unchanged.

Compare the resulting root-source and xylem trajectories.

## Exercise C — leaf-anatomy sensitivity

Change one field in:

```
P.leafAnatomyConditioning
```

and compare the leaf-water and plant-wax responses.


## Exercise D — instantaneous versus standing wax

Run the full reference model, then use the standing-wax integration workflow to compare the instantaneous plant-wax signal with the accumulated standing-wax response.

## Exercise E — time-resolved inverse information

Run:

```
Rinv = test_METCUT_manuscript_inverse_recovery;
```

Compare the final-only and time-resolved uncertainty widths.

---
