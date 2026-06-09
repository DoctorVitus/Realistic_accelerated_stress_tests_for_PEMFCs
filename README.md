# Realistic accelerated stress tests for PEM fuel cells

Dataset and analysis repository for **Realistic accelerated stress tests for PEM fuel cells: Test procedure development based on standardized automotive driving cycles**.

## Repository layout

- `data/`: original dataset files, preserved without editing.
- `code/main.py`: reproducible analysis pipeline.
- `results/`: generated tables and figures.

## Dataset description

Project name: RealAgeing

Type of measurement:

- `PC`: polarization curve
- `EIS`: electrochemical impedance spectroscopy

Reference numbers specify operating conditions. See `data/.../OperatingParameters.xlsx` for details.

For polarization curves:

- `BreakIn`: directly after break-in; no cycles driven yet.
- `xxCycles`: after `xx` cycles.

For EIS:

- Current density is encoded in the filename, for example `0_3` and `0_7` for 0.3 and 0.7 A/cm^2.

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

## Citation

This is the dataset to the published article "Realistic accelerated stress tests for PEM fuel cells: Test procedure development based on standardized automotive driving cycles" (DOI: 10.1016/j.ijhydene.2023.08.292), in which the degradation of two commercial PEM fuel cell stacks was analyzed.

Please cite:

P. Thiele, Y. Yang, S. Dirkes, M. Wick, S. Pischinger, Realistic accelerated stress tests for PEM fuel cells: Test procedure development based on standardized automotive driving cycles, International Journal of Hydrogen Energy 52 (Part D) (2024) 1065-1080, https://doi.org/10.1016/j.ijhydene.2023.08.292.

## Source

- https://zenodo.org/records/13166135

