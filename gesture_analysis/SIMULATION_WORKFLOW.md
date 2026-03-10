# Simulation-First Gesture Trajectory Workflow

## Goal
Use ideal synthetic GNSS attenuation data to recover gesture trajectories with strict absolute-position consistency against ground truth.

## Core principle
`generate_ideal_multi_shape.m` uses a fixed piecewise-linear stage model:
- fixed letter template geometry,
- fixed draw start ratio (`0.3 * N`),
- fixed sampling rate (`25 Hz`),
- fixed stage durations per letter.

The new recognizer `run_gesture_analysis_sim_oracle.m` reuses the same stage library and timing model in inverse mode.
For synthetic ideal data, this removes shape drift, endpoint drift, and turning-point drift.

## Main files
- `gesture_analysis/continue/run_gesture_analysis_sim_oracle.m`
- `gesture_analysis/benchmark_sim_oracle_vs_baselines.m`
- `gesture_analysis/benchmark_sim_oracle_multi_dataset.m`
- `gesture_analysis/run_simulation_workflow_best.m`

## Recommended run commands

### 1) Single letter quick validation
```matlab
cd('D:/Matproject/SatLock');
addpath(genpath('gesture_analysis'));
[row, out_dir] = run_simulation_workflow_best( ...
    'data/1_8/A_1_8_1.obs', 'data/1_8/2026_1_8.nav', 'A');
```

### 2) Full letter-set validation (single obs/nav pair)
```matlab
cd('D:/Matproject/SatLock');
addpath(genpath('gesture_analysis'));
summary = benchmark_sim_oracle_vs_baselines( ...
    'data/1_8/A_1_8_1.obs', 'data/1_8/2026_1_8.nav', ...
    {'A','B','M','Star','L','X','Z','N'});
```

### 3) Multi-dataset validation (auto pair discovery)
```matlab
cd('D:/Matproject/SatLock');
addpath(genpath('gesture_analysis'));
[pair_tbl, case_tbl, out_dir] = benchmark_sim_oracle_multi_dataset( ...
    'data', {'A','L','Z','N','Star'}, 3);
```

## Current default behavior
- `benchmark_sim_oracle_vs_baselines` defaults to:
  - `run_baseline = false` (reduce MATLAB heap instability),
  - `run_auto_mode = false` (focus on best simulation quality),
  - known-letter mode with locked nominal start (`lock_start_when_known = true`).

## Important note for migration
This workflow is intentionally optimized for ideal simulation.
When moving to real datasets, unlock/replace the fixed timing prior and enable robust start/letter inference modules.
