# system_polt

This folder is a standalone workflow for the principle-oriented schematic figures in
`III. SYSTEM DESIGN` of `StarDial.pdf`.

It is intentionally separated from:

- `gesture_analysis/gesture_test.m`
- `gesture_analysis/scientific_graphing/export_paper_figures_data_driven.m`

Its purpose is different from the main evaluation pipeline:

- the main pipeline exports result figures for experiments and metrics
- `system_polt` exports explanatory schematic figures for system principles

## Current Entry Point

Run in MATLAB:

```matlab
[manifest_tbl, out_dir, demo] = generate_system_design_figures();
```

The current generator draws theory-oriented schematic illustrations instead of
reading one real observation day.

## Current Outputs

The script currently exports three principle figures:

- `A_signal_preprocessing`
- `B_diffraction_feature_extraction`
- `C_gesture_plane_geometric_inversion`

Outputs are written to:

```text
system_polt/results/system_design_figures_<timestamp>/
  png/
  pdf/
  fig/
  system_design_manifest.csv
```

## Notes

- The exported figures are intentionally idealized and are meant for paper
  explanation rather than quantitative evidence.
- The B-section feature plot illustrates how one continuous diffraction event
  is converted into amplitude, oscillation, continuity, and fused saliency.
- The C-section figure illustrates how active-contact regions and silent
  exclusions constrain a feasible gesture trajectory.
- This folder is designed so later principle figures can be added without
  touching the main paper-figure exporter.
