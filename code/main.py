from __future__ import annotations

import math
import re
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


def ensure_dirs() -> None:
    for path in [
        TABLE_DIR,
        FIG_DIR / "eis_nyquist",
        FIG_DIR / "eis_trends",
        FIG_DIR / "polarization_curves",
    ]:
        path.mkdir(parents=True, exist_ok=True)


def dataset_dir() -> Path:
    candidates = [p for p in DATA_ROOT.iterdir() if p.is_dir()]
    if not candidates:
        raise FileNotFoundError(f"No dataset directory found under {DATA_ROOT}")
    return candidates[0]


def parse_cycles(label: str) -> int:
    if "breakin" in label.lower() or label.startswith("0_"):
        return 0
    match = re.search(r"(\d+)\s*Cycles", label, flags=re.IGNORECASE)
    return int(match.group(1)) if match else 0


def parse_eis_filename(path: Path) -> tuple[str, float]:
    stem = path.stem
    match = re.search(r"EIS-([A-Za-z0-9]+)-(\d+)_(\d+)$", stem)
    if not match:
        return "unknown", math.nan
    condition = match.group(1)
    current_density = float(f"{match.group(2)}.{match.group(3)}")
    return condition, current_density


def circle_fit(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(x) < 5:
        return {"ok": 0.0}

    a = np.column_stack([x, y, np.ones_like(x)])
    b = -(x**2 + y**2)
    try:
        coef, *_ = np.linalg.lstsq(a, b, rcond=None)
    except np.linalg.LinAlgError:
        return {"ok": 0.0}

    aa, bb, cc = coef
    cx = -aa / 2.0
    cy = -bb / 2.0
    radius_sq = cx**2 + cy**2 - cc
    if radius_sq <= 0:
        return {"ok": 0.0}

    radius = math.sqrt(radius_sq)
    root_sq = radius**2 - cy**2
    if root_sq > 0:
        span = 2.0 * math.sqrt(root_sq)
        left = cx - math.sqrt(root_sq)
        right = cx + math.sqrt(root_sq)
    else:
        span = 2.0 * radius
        left = cx - radius
        right = cx + radius

    residual = np.sqrt((x - cx) ** 2 + (y - cy) ** 2) - radius
    rmse = float(np.sqrt(np.mean(residual**2)))
    return {
        "ok": 1.0,
        "center_x": float(cx),
        "center_y": float(cy),
        "radius": float(radius),
        "x_left": float(left),
        "x_right": float(right),
        "x_span": float(abs(span)),
        "rmse": rmse,
    }


def robust_eis_arrays(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    clean = df[["Z_Re / Ohm", "Z_Im / Ohm"]].apply(pd.to_numeric, errors="coerce").dropna()
    x = clean["Z_Re / Ohm"].to_numpy(dtype=float)
    y = -clean["Z_Im / Ohm"].to_numpy(dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(x) >= 20:
        x_q1, x_q3 = np.quantile(x, [0.25, 0.75])
        y_q1, y_q3 = np.quantile(y, [0.25, 0.75])
        x_iqr = max(x_q3 - x_q1, np.finfo(float).eps)
        y_iqr = max(y_q3 - y_q1, np.finfo(float).eps)
        robust_mask = (
            (x >= x_q1 - 3.0 * x_iqr)
            & (x <= x_q3 + 3.0 * x_iqr)
            & (y >= y_q1 - 3.0 * y_iqr)
            & (y <= y_q3 + 3.0 * y_iqr)
        )
        x = x[robust_mask]
        y = y[robust_mask]
    return x, y


def split_two_arcs(x: np.ndarray, y: np.ndarray) -> int:
    order = np.argsort(x)
    xs = x[order]
    ys = y[order]
    if len(xs) < 12:
        return max(3, len(xs) // 2)

    kernel = np.ones(5) / 5.0
    smooth = np.convolve(ys, kernel, mode="same")
    ymax = float(np.nanmax(smooth))
    candidates: list[int] = []
    for idx in range(1, len(smooth) - 1):
        if smooth[idx] >= smooth[idx - 1] and smooth[idx] >= smooth[idx + 1] and smooth[idx] > ymax * 0.04:
            candidates.append(idx)

    if len(candidates) >= 2:
        first_peak = candidates[0]
        second_peak = max(candidates[1:], key=lambda i: smooth[i])
        if first_peak > second_peak:
            first_peak, second_peak = second_peak, first_peak
        valley = first_peak + int(np.argmin(smooth[first_peak : second_peak + 1]))
        return int(order[valley])

    return int(order[max(4, min(len(xs) - 5, int(len(xs) * 0.22)))])


def fit_eis_sheet(df: pd.DataFrame) -> tuple[dict[str, float | str], pd.DataFrame]:
    x_col = "Z_Re / Ohm"
    x, y = robust_eis_arrays(df)
    if len(x) < 12:
        return {"fit_status": "too_few_points"}, pd.DataFrame()

    split_original = split_two_arcs(x, y)
    order = np.argsort(x)
    split_pos = int(np.where(order == split_original)[0][0])
    split_pos = max(4, min(len(order) - 5, split_pos))
    first_idx = order[: split_pos + 1]
    second_idx = order[split_pos:]

    first = circle_fit(x[first_idx], y[first_idx])
    second = circle_fit(x[second_idx], y[second_idx])
    if not first.get("ok") or not second.get("ok"):
        return {"fit_status": "fit_failed"}, pd.DataFrame()

    params = {
        "fit_status": "ok",
        "split_z_re_ohm": float(x[split_original]),
        "R_ohm_ohm": max(0.0, float(first["x_left"])),
        "R_anode_ct_ohm": float(first["x_span"]),
        "R_cathode_ct_ohm": float(second["x_span"]),
        "R_total_fit_ohm": float(first["x_span"] + second["x_span"]),
        "anode_fit_rmse": float(first["rmse"]),
        "cathode_fit_rmse": float(second["rmse"]),
        "anode_center_x": float(first["center_x"]),
        "anode_center_y": float(first["center_y"]),
        "anode_radius": float(first["radius"]),
        "cathode_center_x": float(second["center_x"]),
        "cathode_center_y": float(second["center_y"]),
        "cathode_radius": float(second["radius"]),
    }

    def fitted_branch(fit: dict[str, float], idx: np.ndarray, arc_name: str, n: int) -> pd.DataFrame:
        lo = max(float(np.min(x[idx])), float(fit["center_x"] - fit["radius"]))
        hi = min(float(np.max(x[idx])), float(fit["center_x"] + fit["radius"]))
        x_line = np.linspace(lo, hi, n)
        inside = np.maximum(0.0, float(fit["radius"]) ** 2 - (x_line - float(fit["center_x"])) ** 2)
        y_line = float(fit["center_y"]) + np.sqrt(inside)
        return pd.DataFrame({"arc": arc_name, x_col: x_line, "-Z_Im / Ohm": y_line})

    fit_points = pd.concat(
        [
            fitted_branch(first, first_idx, "anode_high_frequency", 140),
            fitted_branch(second, second_idx, "cathode_low_frequency", 180),
        ],
        ignore_index=True,
    )
    return params, fit_points


def analyze_eis(base: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    nyquist_examples: dict[tuple[str, float], list[tuple[Path, str, pd.DataFrame, pd.DataFrame, dict[str, object]]]] = {}

    for path in sorted(base.glob("Stack*/EIS/*/*.xlsx")):
        stack = path.parts[-4]
        cycle_label = path.parts[-2]
        cycles = parse_cycles(cycle_label)
        condition, current_density = parse_eis_filename(path)
        excel = pd.ExcelFile(path)
        for sheet in excel.sheet_names:
            df = pd.read_excel(path, sheet_name=sheet)
            if {"Z_Re / Ohm", "Z_Im / Ohm"}.issubset(df.columns):
                params, fit_points = fit_eis_sheet(df)
                row = {
                    "stack": stack,
                    "cycle_label": cycle_label,
                    "cycles": cycles,
                    "condition": condition,
                    "current_density_A_cm2": current_density,
                    "file": str(path.relative_to(ROOT)),
                    "sheet": sheet,
                    "sheet_label": SHEET_LABELS.get(sheet, sheet),
                }
                row.update(params)
                rows.append(row)
                if sheet == "Sheet1" and params.get("fit_status") == "ok":
                    nyquist_examples.setdefault((stack, current_density), []).append((path, cycle_label, df, fit_points, row))

    eis = pd.DataFrame(rows)
    if eis.empty:
        return eis

    eis.to_csv(TABLE_DIR / "eis_fit_parameters.csv", index=False, encoding="utf-8-sig")
    eis.to_excel(TABLE_DIR / "eis_fit_parameters.xlsx", index=False)

    plot_eis_trends(eis)
    plot_nyquist_examples(nyquist_examples)
    return eis


def plot_eis_trends(eis: pd.DataFrame) -> None:
    ok = eis[eis["fit_status"] == "ok"].copy()
    components = [
        ("R_ohm_ohm", "Ohmic resistance / Ohm"),
        ("R_anode_ct_ohm", "Anode charge-transfer resistance / Ohm"),
        ("R_cathode_ct_ohm", "Cathode charge-transfer resistance / Ohm"),
    ]
    for (stack, current, sheet, sheet_label), group in ok.groupby(
        ["stack", "current_density_A_cm2", "sheet", "sheet_label"], dropna=False
    ):
        fig, axes = plt.subplots(1, 3, figsize=(12, 3.4), constrained_layout=True)
        for ax, (col, ylabel) in zip(axes, components):
            ax.scatter(group["cycles"], group[col], s=18, alpha=0.75)
            mean_line = group.groupby("cycles", as_index=False)[col].mean().sort_values("cycles")
            ax.plot(mean_line["cycles"], mean_line[col], linestyle=":", linewidth=1.4)
            ax.set_xlabel("Number of cycles")
            ax.set_ylabel(ylabel)
            ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.55)
        fig.suptitle(f"{stack} {sheet_label}, {current:g} A/cm^2")
        safe_label = str(sheet_label).replace(" ", "_")
        fig.savefig(FIG_DIR / "eis_trends" / f"{stack}_{safe_label}_{current:g}Acm2_resistance_trends.png")
        plt.close(fig)


def plot_nyquist_examples(examples: dict[tuple[str, float], list[tuple[Path, str, pd.DataFrame, pd.DataFrame, dict[str, object]]]]) -> None:
    for (stack, current), entries in examples.items():
        fig, ax = plt.subplots(figsize=(5.6, 4.4), constrained_layout=True)
        selected = sorted(entries, key=lambda item: parse_cycles(item[1]))
        for path, cycle_label, df, fit_points, _row in selected:
            cycles = parse_cycles(cycle_label)
            condition = _row.get("condition", "")
            label = f"{cycles} cycles, #{condition}"
            x, y = robust_eis_arrays(df)
            ax.scatter(x, y, s=9, alpha=0.65, label=label)
            ax.plot(fit_points["Z_Re / Ohm"], fit_points["-Z_Im / Ohm"], linestyle=":", linewidth=1.0)
        ax.set_xlabel("Z_Re / Ohm")
        ax.set_ylabel("-Z_Im / Ohm")
        ax.set_title(f"{stack} whole stack EIS, {current:g} A/cm^2")
        ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.55)
        ax.legend(ncol=2)
        fig.savefig(FIG_DIR / "eis_nyquist" / f"{stack}_{current:g}Acm2_whole_stack_nyquist_fit.png")
        plt.close(fig)


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


def analyze_polcurves(base: Path) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for path in sorted(base.glob("Stack*/PolCurves/*.xlsx")):
        stack = path.parts[-3]
        match = re.search(r"PC_(.+?)_1$", path.stem)
        cycle_label = match.group(1) if match else path.stem
        cycles = parse_cycles(cycle_label)
        df = read_polcurve(path)
        df.insert(0, "file", str(path.relative_to(ROOT)))
        df.insert(0, "cycles", cycles)
        df.insert(0, "cycle_label", cycle_label)
        df.insert(0, "stack", stack)
        rows.append(df)

    pol = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if pol.empty:
        return pol

    pol.to_csv(TABLE_DIR / "polarization_curve_points.csv", index=False, encoding="utf-8-sig")
    pol.to_excel(TABLE_DIR / "polarization_curve_points.xlsx", index=False)
    plot_polcurves(pol)
    return pol


def plot_polcurves(pol: pd.DataFrame) -> None:
    x_col = "Cmon_Current (A/cm^2)"
    y_candidates = ["Average_Voltage (mV)", "Minimum_Voltage (mV)", "Maximum_Voltage (mV)"]
    for stack, group in pol.groupby("stack"):
        fig, ax = plt.subplots(figsize=(5.8, 4.2), constrained_layout=True)
        for cycle_label, curve in group.groupby("cycle_label"):
            curve = curve.sort_values(x_col)
            ax.scatter(curve[x_col], curve["Average_Voltage (mV)"], s=16, label=str(cycle_label), alpha=0.8)
            ax.plot(curve[x_col], curve["Average_Voltage (mV)"], linestyle=":", linewidth=1.2)
        ax.set_xlabel(x_col)
        ax.set_ylabel("Average_Voltage (mV)")
        ax.set_title(f"{stack} polarization curves")
        ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.55)
        ax.legend()
        fig.savefig(FIG_DIR / "polarization_curves" / f"{stack}_average_voltage_polarization_curves.png")
        plt.close(fig)

        for y_col in y_candidates:
            if y_col not in group.columns:
                continue
            summary = group.groupby("cycles", as_index=False)[y_col].mean().sort_values("cycles")
            fig, ax = plt.subplots(figsize=(5.2, 3.7), constrained_layout=True)
            ax.scatter(summary["cycles"], summary[y_col], s=22)
            ax.plot(summary["cycles"], summary[y_col], linestyle=":", linewidth=1.3)
            ax.set_xlabel("Number of cycles")
            ax.set_ylabel(y_col)
            ax.set_title(f"{stack} mean {y_col}")
            ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.55)
            safe_y = y_col.replace(" ", "_").replace("/", "_").replace("^", "")
            fig.savefig(FIG_DIR / "polarization_curves" / f"{stack}_{safe_y}_trend.png")
            plt.close(fig)


def write_summary(eis: pd.DataFrame, pol: pd.DataFrame) -> None:
    lines = [
        "# PEMFC Ageing Dataset Analysis Results",
        "",
        "This folder contains generated tables and figures from the original PEM fuel-cell ageing dataset.",
        "",
        "## EIS fitting convention",
        "",
        "- Nyquist plots use `Z_Re / Ohm` versus `-Z_Im / Ohm`.",
        "- `R_ohm_ohm` is estimated from the high-frequency x-axis intercept of the first fitted semicircle.",
        "- `R_anode_ct_ohm` is assigned to the high-frequency semicircle diameter.",
        "- `R_cathode_ct_ohm` is assigned to the low-frequency semicircle diameter.",
        "- Fits use algebraic circle fitting, so the script does not require SciPy.",
        "",
        "## Generated tables",
        "",
        "- `tables/eis_fit_parameters.csv` and `.xlsx`",
        "- `tables/polarization_curve_points.csv` and `.xlsx`",
        "",
        "## Generated figures",
        "",
        "- `figures/eis_nyquist`: whole-stack Nyquist scatter plots with dotted fitted semicircles.",
        "- `figures/eis_trends`: resistance trends by stack, current density, and sheet.",
        "- `figures/polarization_curves`: polarization curves and voltage trends.",
        "",
        f"EIS fit rows: {len(eis)}",
        f"Polarization curve rows: {len(pol)}",
    ]
    (RESULTS_ROOT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()
    base = dataset_dir()
    eis = analyze_eis(base)
    pol = analyze_polcurves(base)
    write_summary(eis, pol)
    print(f"EIS rows: {len(eis)}")
    print(f"Polarization rows: {len(pol)}")
    print(f"Results written to: {RESULTS_ROOT}")


if __name__ == "__main__":
    main()
