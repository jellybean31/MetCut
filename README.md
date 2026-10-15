# MetCut

**MetCut** is a process-oriented model of hydrogen-isotope propagation from environmental water to plant wax, implemented in MATLAB. It combines soil-water and root uptake processes, plant hydraulics, xylem isotope transport and storage, leaf-water evaporative enrichment, intracellular water exchange, biochemical hydrogen incorporation, and standing plant-wax accumulation.

The model tracks δD through intermediate water and biochemical states and predicts instantaneous and accumulated plant-wax δD for the C27, C29, C31, and C33 alkane homologues. The package also includes sensitivity, tracer, and synthetic inverse-recovery examples.

## Citation
If you use MetCut for future publications, please cite: [`pending...`]

## Add MetCut to the MATLAB Path
To use MetCut functions and examples from any MATLAB working directory, add the functions/ and examples/ folders to your MATLAB search path.

Open MATLAB and run:

```matlab
modelRoot = '/path/to/MetCut/folder';

addpath(fullfile(modelRoot, 'functions'));
addpath(fullfile(modelRoot, 'examples'));

savepath;
```

Replace /path/to/MetCut/folder with the path to your downloaded or cloned MetCut directory.
Saving the path: The savepath command keeps these changes across MATLAB sessions, so you do not need to repeat the setup each time MATLAB starts. Alternatively, use Home → Set Path → Save in the MATLAB interface.

To verify installation, run:

```matlab
which run_full_model_from_hydraulic_checkpoint
```

MATLAB should return the full path to the corresponding function in MetCut/functions/.
Working directory: Once added, MetCut functions can be called from any MATLAB working directory. However, some example and validation scripts use relative paths to locate data or checkpoints. These scripts should be executed from the MetCut root directory unless otherwise specified in their instructions.
Note: If you move the MetCut directory after installation, update the saved MATLAB path to reflect its new location.

## Repository contents

```text
MetCut/
├── functions/          Model functions and numerical utilities
├── data/
│   ├── reference/      Reference model checkpoints and validation data
│   └── DATA_SOURCES.md Source datasets and literature used in development
├── examples/
│   ├── METCUT_v1_Tutorial.md
│   ├── example_01_reference_forward_run.m
│   ├── example_02_parameterized_sensitivity.m
│   ├── example_03_precipitation_tracer.m
│   ├── example_04_standing_wax_integration.m
│   ├── example_05_inverse_recovery.m
│   └── run_all_examples.m
└── tests/
    ├── unit/           Numerical and component tests
    ├── integration/    Coupled-model and standing-wax tests
    └── manuscript/     Tests associated with reported model examples
```

## Getting started

MATLAB is required. Start MATLAB in the root `MetCut` directory and run:

```matlab
root = pwd;
addpath(fullfile(root, 'functions'), '-begin');
addpath(fullfile(root, 'examples'), '-begin');

R = example_01_reference_forward_run;
```

This example loads the supplied eight-day hydraulic reference checkpoint, propagates the isotope calculations, and plots intermediate water and instantaneous plant-wax δD trajectories. **It reuses the saved hydraulic solution.**

To run the reference forward model directly:

```matlab
hydPath = fullfile(root, 'data', 'reference', ...
    'hydraulic_checkpoint_76cell_8day_hourly_corrected_mass_conservative_storagefixed.mat');

[M, B] = run_full_model_from_hydraulic_checkpoint(hydPath);
```

`M` contains the modeled states and trajectories; `B` contains the isotope-fractionation budget and related diagnostics.

## Examples and tests

The five examples cover the reference forward calculation, a parameterized precipitation-isotope sensitivity test, a precipitation tracer, standing-wax integration, and a synthetic inverse-recovery demonstration. Run them individually using the corresponding `example_0*.m` functions, or run `run_all_examples` from the `examples/` folder on the MATLAB path.

For validation, add the appropriate test directory to the MATLAB path and run one of the included test suites:

```matlab
addpath(fullfile(root, 'tests', 'manuscript'), '-begin');
R = run_METCUT_manuscript_tests;
```

Some tests and examples, especially coupled-hydraulic and inverse-recovery calculations, can be computationally expensive. The supplied reference checkpoints allow many demonstrations to run without repeating the hydraulic solution.

## Documentation and data

See [`examples/METCUT_v1_Tutorial.md`] for module descriptions, additional MATLAB commands, and worked examples.

The `data/reference/` directory contains the reference files required by the packaged examples. [`data/DATA_SOURCES.md`] documents published datasets and other sources used during model development.
