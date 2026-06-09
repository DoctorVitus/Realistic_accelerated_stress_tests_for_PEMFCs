# Realistic Accelerated Stress Tests for PEMFCs

This repository organizes the original PEM fuel-cell ageing dataset and generated analysis outputs.

## Repository layout

- `data/`: original dataset files, preserved without editing.
- `code/main.py`: reproducible analysis pipeline.
- `results/`: generated tables and figures.

## Analysis summary

The analysis script reads EIS and polarization-curve Excel files, generates scatter plots with dotted guide/fit lines, and writes tabular summaries.

For EIS, Nyquist data are plotted as `Z_Re / Ohm` versus `-Z_Im / Ohm`. The two visible semicircles are estimated with algebraic circle fitting:

- `R_ohm_ohm`: high-frequency x-axis intercept estimate.
- `R_anode_ct_ohm`: high-frequency semicircle diameter.
- `R_cathode_ct_ohm`: low-frequency semicircle diameter.

All generated figures use Times New Roman font settings in Matplotlib.

## Reproduce results

```bash
python code/main.py
```

