from __future__ import annotations

import re
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT / "data"
RESULTS_ROOT = ROOT / "results"
TABLE_DIR = RESULTS_ROOT / "tables"
FIG_DIR = RESULTS_ROOT / "figures"

SHEET_LABELS = {
    "Sheet1": "Whole Stack",
    "Sheet2": "Cell 01",
    "Sheet3": "Cell 03",
    "Sheet4": "Cell 05",
    "Sheet5": "Cell 10",
    "Sheet6": "Cell 15",
    "Sheet7": "Cell 20",
    "Sheet8": "Cell 22",
    "Sheet9": "Cell 24",
}

plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "legend.fontsize": 8,
        "figure.dpi": 150,
        "savefig.dpi": 300,
    }
)


def reset_results() -> None:
    if RESULTS_ROOT.exists():
        shutil.rmtree(RESULTS_ROOT, ignore_errors=True)
    for path in [
        TABLE_DIR,
        FIG_DIR / "eis" / "individual",
        FIG_DIR / "eis" / "by_condition",
        FIG_DIR / "polarization" / "individual",
        FIG_DIR / "polarization" / "by_stack",
    ]:
        path.mkdir(parents=True, exist_ok=True)


def parse_cycles(label: str) -> int:
    if "breakin" in label.lower() or label.startswith("0_"):
        return 0
    match = re.search(r"(\d+)\s*Cycles", label, flags=re.IGNORECASE)
    return int(match.group(1)) if match else 0


def parse_eis_filename(path: Path) -> tuple[str, float]:
    match = re.search(r"EIS-([A-Za-z0-9]+)-(\d+)_(\d+)$", path.stem)
    if not match:
        return "unknown", np.nan
    return match.group(1), float(f"{match.group(2)}.{match.group(3)}")


def safe_name(value: object) -> str:
    text = str(value)
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    return text.strip("_")


def robust_eis_data(df: pd.DataFrame) -> pd.DataFrame:
    required = ["Frequency / Hz", "Z_Re / Ohm", "Z_Im / Ohm"]
    clean = df[required].apply(pd.to_numeric, errors="coerce").dropna()
    clean = clean.replace([np.inf, -np.inf], np.nan).dropna()
    if len(clean) >= 20:
        x = clean["Z_Re / Ohm"]
        y = -clean["Z_Im / Ohm"]
        x_q1, x_q3 = np.quantile(x, [0.25, 0.75])
        y_q1, y_q3 = np.quantile(y, [0.25, 0.75])
        x_iqr = max(float(x_q3 - x_q1), np.finfo(float).eps)
        y_iqr = max(float(y_q3 - y_q1), np.finfo(float).eps)
        mask = (
            (x >= x_q1 - 3.0 * x_iqr)
            & (x <= x_q3 + 3.0 * x_iqr)
            & (y >= y_q1 - 3.0 * y_iqr)
            & (y <= y_q3 + 3.0 * y_iqr)
        )
        clean = clean[mask].copy()
    clean["-Z_Im / Ohm"] = -clean["Z_Im / Ohm"]
    return clean.sort_values("Frequency / Hz", ascending=False)


def scatter_with_guide(
    ax: plt.Axes,
    x: pd.Series,
    y: pd.Series,
    *,
    size: float,
    alpha: float,
    line_width: float,
    label: str | None = None,
) -> None:
    line = ax.plot(x, y, linestyle=":", linewidth=line_width, alpha=alpha)[0]
    ax.scatter(
        x,
        y,
        s=size,
        facecolors="none",
        edgecolors=line.get_color(),
        linewidths=0.8,
        alpha=alpha,
        label=label,
    )


def plot_eis_single(df: pd.DataFrame, meta: dict[str, object]) -> Path:
    fig, ax = plt.subplots(figsize=(5.0, 3.8), constrained_layout=True)
    scatter_with_guide(
        ax,
        df["Z_Re / Ohm"],
        df["-Z_Im / Ohm"],
        size=15,
        alpha=0.78,
        line_width=1.2,
    )
    ax.set_xlabel("Z_Re / Ohm")
    ax.set_ylabel("-Z_Im / Ohm")
    ax.set_title(
        f"{meta['stack']} {meta['sheet_label']}, {meta['cycles']} cycles, "
        f"#{meta['condition']}, {meta['current_density_A_cm2']:g} A/cm^2"
    )
    ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.55)

    out_dir = (
        FIG_DIR
        / "eis"
        / "individual"
        / safe_name(meta["stack"])
        / safe_name(meta["sheet_label"])
        / f"{meta['current_density_A_cm2']:g}Acm2"
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / (
        f"{safe_name(meta['cycle_label'])}_condition_{safe_name(meta['condition'])}_"
        f"{safe_name(Path(str(meta['file'])).stem)}.png"
    )
    fig.savefig(out)
    plt.close(fig)
    return out


def plot_eis_group(items: list[tuple[pd.DataFrame, dict[str, object]]], key: tuple[object, ...]) -> Path:
    stack, sheet_label, condition, current = key
    fig, ax = plt.subplots(figsize=(5.8, 4.2), constrained_layout=True)
    for df, meta in sorted(items, key=lambda item: (int(item[1]["cycles"]), str(item[1]["file"]))):
        label = f"{int(meta['cycles'])} cycles"
        scatter_with_guide(
            ax,
            df["Z_Re / Ohm"],
            df["-Z_Im / Ohm"],
            size=10,
            alpha=0.65,
            line_width=1.0,
            label=label,
        )
    ax.set_xlabel("Z_Re / Ohm")
    ax.set_ylabel("-Z_Im / Ohm")
    ax.set_title(f"{stack} {sheet_label}, condition {condition}, {current:g} A/cm^2")
    ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.55)
    ax.legend(ncol=2)

    out_dir = FIG_DIR / "eis" / "by_condition" / safe_name(stack) / safe_name(sheet_label)
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"condition_{safe_name(condition)}_{current:g}Acm2.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def analyze_eis() -> pd.DataFrame:
    records: list[dict[str, object]] = []
    grouped_data: dict[tuple[object, ...], list[tuple[pd.DataFrame, dict[str, object]]]] = {}

    for path in sorted(DATA_ROOT.glob("Stack*/EIS/*/*.xlsx")):
        stack = path.parts[-4]
        cycle_label = path.parts[-2]
        cycles = parse_cycles(cycle_label)
        condition, current = parse_eis_filename(path)
        excel = pd.ExcelFile(path)
        for sheet in excel.sheet_names:
            raw = pd.read_excel(path, sheet_name=sheet)
            if not {"Frequency / Hz", "Z_Re / Ohm", "Z_Im / Ohm"}.issubset(raw.columns):
                continue
            df = robust_eis_data(raw)
            sheet_label = SHEET_LABELS.get(sheet, sheet)
            meta = {
                "stack": stack,
                "cycle_label": cycle_label,
                "cycles": cycles,
                "condition": condition,
                "current_density_A_cm2": current,
                "file": str(path.relative_to(ROOT)),
                "sheet": sheet,
                "sheet_label": sheet_label,
            }
            fig_path = plot_eis_single(df, meta)
            meta["figure"] = str(fig_path.relative_to(ROOT))
            records.append(meta)
            group_key = (stack, sheet_label, condition, current)
            grouped_data.setdefault(group_key, []).append((df, meta))

    index = pd.DataFrame(records)
    if index.empty:
        return index

    group_records: list[dict[str, object]] = []
    keys = ["stack", "sheet_label", "condition", "current_density_A_cm2"]
    for key, items in grouped_data.items():
        fig_path = plot_eis_group(items, key)
        group_records.append(
            {
                "stack": key[0],
                "sheet_label": key[1],
                "condition": key[2],
                "current_density_A_cm2": key[3],
                "figure": str(fig_path.relative_to(ROOT)),
                "curves": len(items),
            }
        )

    groups = pd.DataFrame(group_records).sort_values(keys)
    index.to_csv(TABLE_DIR / "eis_plot_index.csv", index=False, encoding="utf-8-sig")
    groups.to_csv(TABLE_DIR / "eis_condition_plot_index.csv", index=False, encoding="utf-8-sig")
    return index


def read_polcurve(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path, header=None)
    header = raw.iloc[1].astype(str).tolist()
    units = raw.iloc[2].astype(str).tolist()
    data = raw.iloc[3:].copy()
    data.columns = header
    data = data.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    for col, unit in zip(header, units):
        if unit != "nan":
            data.rename(columns={col: f"{col} ({unit})"}, inplace=True)
    return data


def plot_pol_single(df: pd.DataFrame, meta: dict[str, object]) -> Path:
    x_col = "Cmon_Current (A/cm^2)"
    voltage_cols = [col for col in df.columns if col.endswith("(mV)") and col != x_col]
    fig, ax = plt.subplots(figsize=(5.4, 4.0), constrained_layout=True)
    for col in voltage_cols:
        if col in {"Minimum_Voltage (mV)", "Average_Voltage (mV)", "Maximum_Voltage (mV)"}:
            line_width = 1.3
            alpha = 0.9
        else:
            line_width = 0.9
            alpha = 0.55
        curve = df[[x_col, col]].dropna().sort_values(x_col)
        scatter_with_guide(
            ax,
            curve[x_col],
            curve[col],
            size=12,
            alpha=alpha,
            line_width=line_width,
            label=col.replace(" (mV)", ""),
        )
    ax.set_xlabel(x_col)
    ax.set_ylabel("Voltage (mV)")
    ax.set_title(f"{meta['stack']} polarization, {meta['cycle_label']}")
    ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.55)
    ax.legend(ncol=2)

    out_dir = FIG_DIR / "polarization" / "individual" / safe_name(meta["stack"])
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{safe_name(meta['cycle_label'])}_{safe_name(Path(str(meta['file'])).stem)}.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def plot_pol_stack(group: pd.DataFrame, stack: str) -> Path:
    x_col = "Cmon_Current (A/cm^2)"
    y_col = "Average_Voltage (mV)"
    fig, ax = plt.subplots(figsize=(5.8, 4.2), constrained_layout=True)
    for cycle_label, curve in group.groupby("cycle_label"):
        curve = curve[[x_col, y_col]].dropna().sort_values(x_col)
        scatter_with_guide(
            ax,
            curve[x_col],
            curve[y_col],
            size=16,
            alpha=0.75,
            line_width=1.2,
            label=str(cycle_label),
        )
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title(f"{stack} average-voltage polarization curves")
    ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.55)
    ax.legend()

    out_dir = FIG_DIR / "polarization" / "by_stack"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{safe_name(stack)}_average_voltage.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def analyze_polarization() -> tuple[pd.DataFrame, pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    records: list[dict[str, object]] = []
    for path in sorted(DATA_ROOT.glob("Stack*/PolCurves/*.xlsx")):
        stack = path.parts[-3]
        match = re.search(r"PC_(.+?)_1$", path.stem)
        cycle_label = match.group(1) if match else path.stem
        cycles = parse_cycles(cycle_label)
        df = read_polcurve(path)
        meta = {
            "stack": stack,
            "cycle_label": cycle_label,
            "cycles": cycles,
            "file": str(path.relative_to(ROOT)),
        }
        fig_path = plot_pol_single(df, meta)
        meta["figure"] = str(fig_path.relative_to(ROOT))
        records.append(meta)
        df.insert(0, "file", str(path.relative_to(ROOT)))
        df.insert(0, "cycles", cycles)
        df.insert(0, "cycle_label", cycle_label)
        df.insert(0, "stack", stack)
        frames.append(df)

    points = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    index = pd.DataFrame(records)
    if points.empty:
        return points, index

    stack_records = []
    for stack, group in points.groupby("stack"):
        fig_path = plot_pol_stack(group, stack)
        stack_records.append({"stack": stack, "figure": str(fig_path.relative_to(ROOT)), "curves": group["cycle_label"].nunique()})

    points.to_csv(TABLE_DIR / "polarization_curve_points.csv", index=False, encoding="utf-8-sig")
    index.to_csv(TABLE_DIR / "polarization_plot_index.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(stack_records).to_csv(TABLE_DIR / "polarization_stack_plot_index.csv", index=False, encoding="utf-8-sig")
    return points, index


def operating_conditions_table() -> pd.DataFrame:
    path = DATA_ROOT / "OperatingParameters.xlsx"
    raw = pd.read_excel(path, header=None)
    header = raw.iloc[0].tolist()
    units = raw.iloc[1].tolist()
    data = raw.iloc[2:, :5].copy()
    data.columns = header[:5]
    data = data.dropna(how="all")
    for col, unit in zip(header[:5], units[:5]):
        if isinstance(unit, str) and unit not in {"#", "/"}:
            data.rename(columns={col: f"{col} ({unit})"}, inplace=True)
    data.to_csv(TABLE_DIR / "operating_conditions.csv", index=False, encoding="utf-8-sig")
    return data


def write_results_readme(eis_index: pd.DataFrame, pol_index: pd.DataFrame) -> None:
    text = f"""# Generated Results

This folder contains visualization-only outputs generated from the original dataset.

- EIS individual plots: {len(eis_index)}
- Polarization individual plots: {len(pol_index)}

No resistance extraction or curve-fitting-derived parameter table is included.

## Main folders

- `figures/eis/individual/`: one Nyquist plot per EIS workbook sheet.
- `figures/eis/by_condition/`: grouped Nyquist plots by stack, sheet, condition, and current density.
- `figures/polarization/individual/`: one polarization plot per polarization workbook.
- `figures/polarization/by_stack/`: grouped average-voltage polarization curves by stack.
- `tables/`: plot index tables and cleaned plotting points.
"""
    (RESULTS_ROOT / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    reset_results()
    operating_conditions_table()
    eis_index = analyze_eis()
    _, pol_index = analyze_polarization()
    write_results_readme(eis_index, pol_index)
    print(f"EIS plots: {len(eis_index)}")
    print(f"Polarization plots: {len(pol_index)}")
    print(f"Results written to: {RESULTS_ROOT}")


if __name__ == "__main__":
    main()
