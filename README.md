# Realistic Accelerated Stress Tests for PEM Fuel Cells

This repository contains the RealAgeing PEMFC dataset and a visualization-only workflow for quick inspection of EIS and polarization-curve data.

## Repository Structure

```text
data/     Original dataset files, preserved without editing.
code/     Reproducible Python visualization script.
results/  Generated figures and plot index tables.
```

The data folder is kept clean at the repository root:

```text
data/
  FileNameCode.txt
  OperatingParameters.xlsx
  Stack1/
  Stack2/
```

## Quick Run

```bash
python code/main.py
```

The script regenerates `results/` from the Excel files under `data/`.

## Visualization Notes

- EIS plots are raw Nyquist-style visualizations: `Z_Re / Ohm` vs `-Z_Im / Ohm`.
- No semicircle fitting, resistance extraction, or fitted resistance tracking is included.
- All plots use scatter markers plus dotted guide lines.
- Axis labels follow original column names where applicable.
- Matplotlib is configured to use Times New Roman.

## Operating Condition Numbers

The number in filenames such as `Real_Ageing_EIS-1-0_3.xlsx` is the operating-condition reference number from `data/OperatingParameters.xlsx`.

| Condition | Temperature (C) | Air pressure (bar) | Air RH (%) | Air stoichiometry |
| --- | ---: | ---: | ---: | ---: |
| 1 | 60 | 1.50 | 70 | 2 |
| 2 | 60 | 1.75 | 70 | 2 |
| 6 | 60 | 1.25 | 70 | 2 |
| 7 | 60 | 1.25 | 70 | 3 |

The full condition table is saved at `results/tables/operating_conditions.csv`.

## EIS By Condition

Grouped EIS figures are saved under `results/figures/eis/by_condition/`. Each plot groups curves by cycle count for the same stack, sheet, condition, and current density.

### Condition 1

#### Stack1, Whole Stack, 0.3 A/cm^2

![Stack1 condition 1 0.3 A/cm2](results/figures/eis/by_condition/Stack1/Whole_Stack/condition_1_0.3Acm2.png)

#### Stack1, Whole Stack, 0.7 A/cm^2

![Stack1 condition 1 0.7 A/cm2](results/figures/eis/by_condition/Stack1/Whole_Stack/condition_1_0.7Acm2.png)

#### Stack2, Whole Stack, 0.3 A/cm^2

![Stack2 condition 1 0.3 A/cm2](results/figures/eis/by_condition/Stack2/Whole_Stack/condition_1_0.3Acm2.png)

#### Stack2, Whole Stack, 0.7 A/cm^2

![Stack2 condition 1 0.7 A/cm2](results/figures/eis/by_condition/Stack2/Whole_Stack/condition_1_0.7Acm2.png)

### Condition 2

#### Stack1, Whole Stack, 0.3 A/cm^2

![Stack1 condition 2 0.3 A/cm2](results/figures/eis/by_condition/Stack1/Whole_Stack/condition_2_0.3Acm2.png)

#### Stack1, Whole Stack, 0.7 A/cm^2

![Stack1 condition 2 0.7 A/cm2](results/figures/eis/by_condition/Stack1/Whole_Stack/condition_2_0.7Acm2.png)

#### Stack2, Whole Stack, 0.3 A/cm^2

![Stack2 condition 2 0.3 A/cm2](results/figures/eis/by_condition/Stack2/Whole_Stack/condition_2_0.3Acm2.png)

#### Stack2, Whole Stack, 0.7 A/cm^2

![Stack2 condition 2 0.7 A/cm2](results/figures/eis/by_condition/Stack2/Whole_Stack/condition_2_0.7Acm2.png)

### Condition 6

#### Stack1, Whole Stack, 0.3 A/cm^2

![Stack1 condition 6 0.3 A/cm2](results/figures/eis/by_condition/Stack1/Whole_Stack/condition_6_0.3Acm2.png)

#### Stack1, Whole Stack, 0.7 A/cm^2

![Stack1 condition 6 0.7 A/cm2](results/figures/eis/by_condition/Stack1/Whole_Stack/condition_6_0.7Acm2.png)

#### Stack2, Whole Stack, 0.3 A/cm^2

![Stack2 condition 6 0.3 A/cm2](results/figures/eis/by_condition/Stack2/Whole_Stack/condition_6_0.3Acm2.png)

#### Stack2, Whole Stack, 0.7 A/cm^2

![Stack2 condition 6 0.7 A/cm2](results/figures/eis/by_condition/Stack2/Whole_Stack/condition_6_0.7Acm2.png)

### Condition 7

#### Stack1, Whole Stack, 0.3 A/cm^2

![Stack1 condition 7 0.3 A/cm2](results/figures/eis/by_condition/Stack1/Whole_Stack/condition_7_0.3Acm2.png)

#### Stack1, Whole Stack, 0.7 A/cm^2

![Stack1 condition 7 0.7 A/cm2](results/figures/eis/by_condition/Stack1/Whole_Stack/condition_7_0.7Acm2.png)

#### Stack2, Whole Stack, 0.3 A/cm^2

![Stack2 condition 7 0.3 A/cm2](results/figures/eis/by_condition/Stack2/Whole_Stack/condition_7_0.3Acm2.png)

#### Stack2, Whole Stack, 0.7 A/cm^2

![Stack2 condition 7 0.7 A/cm2](results/figures/eis/by_condition/Stack2/Whole_Stack/condition_7_0.7Acm2.png)

## Polarization Curves

Grouped polarization figures are saved under `results/figures/polarization/by_stack/`.

### Stack1

![Stack1 polarization](results/figures/polarization/by_stack/Stack1_average_voltage.png)

### Stack2

![Stack2 polarization](results/figures/polarization/by_stack/Stack2_average_voltage.png)

## Output Inventory

- `results/figures/eis/individual/`: one EIS plot per workbook sheet.
- `results/figures/eis/by_condition/`: grouped EIS plots by condition.
- `results/figures/polarization/individual/`: one plot per polarization workbook.
- `results/figures/polarization/by_stack/`: grouped polarization plots by stack.
- `results/tables/eis_plot_index.csv`: index of all individual EIS figures.
- `results/tables/eis_condition_plot_index.csv`: index of all grouped EIS condition figures.
- `results/tables/polarization_plot_index.csv`: index of all individual polarization figures.

## Citation

This is the dataset to the published article "Realistic accelerated stress tests for PEM fuel cells: Test procedure development based on standardized automotive driving cycles" (DOI: 10.1016/j.ijhydene.2023.08.292).

Please cite:

P. Thiele, Y. Yang, S. Dirkes, M. Wick, S. Pischinger, Realistic accelerated stress tests for PEM fuel cells: Test procedure development based on standardized automotive driving cycles, International Journal of Hydrogen Energy 52 (Part D) (2024) 1065-1080, https://doi.org/10.1016/j.ijhydene.2023.08.292.

## Source

- https://zenodo.org/records/13166135

