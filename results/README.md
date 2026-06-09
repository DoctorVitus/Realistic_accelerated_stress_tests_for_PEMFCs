# PEMFC Ageing Dataset Analysis Results

This folder contains generated tables and figures from the original PEM fuel-cell ageing dataset.

## EIS fitting convention

- Nyquist plots use `Z_Re / Ohm` versus `-Z_Im / Ohm`.
- `R_ohm_ohm` is estimated from the high-frequency x-axis intercept of the first fitted semicircle.
- `R_anode_ct_ohm` is assigned to the high-frequency semicircle diameter.
- `R_cathode_ct_ohm` is assigned to the low-frequency semicircle diameter.
- Fits use algebraic circle fitting, so the script does not require SciPy.

## Generated tables

- `tables/eis_fit_parameters.csv` and `.xlsx`
- `tables/polarization_curve_points.csv` and `.xlsx`

## Generated figures

- `figures/eis_fits`: individual Nyquist scatter plots with dotted two-circle fits.
- `figures/eis_trends`: resistance trends by stack, current density, and sheet.
- `figures/polarization_curves`: polarization curves and voltage trends.

EIS fit rows: 801
Polarization curve rows: 287
