"""
Generate publication-quality figures from Stoixeion fast exports.
Plots:
  1. Within-day Core Jaccard Recurrence Heatmaps (OF1 -> SAMPLE -> TEST -> OF2)
  2. Global Factor Phase Activity Dynamics (event rates across phases)
  3. Consolidated 3x3 Panel: HabL vs SD VEH vs SD CNO across R004, R005, R006.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Set elegant publication aesthetics
plt.rcParams.update({
    "font.sans-serif": "DejaVu Sans",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "figure.titlesize": 14,
    "figure.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})

PHASES_TASK = ["OF1", "SAMPLE", "TEST", "OF2"]
PHASES_HABL = ["OF1", "OF2", "OF3", "OF4"]


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment."""
    p_values = np.asarray(p_values, dtype=float)
    n = len(p_values)
    if n == 0:
        return np.array([])
    order = np.argsort(p_values)
    adjusted = np.empty(n, dtype=float)
    running = 1.0
    for reverse_rank, idx in enumerate(order[::-1], start=1):
        rank = n - reverse_rank + 1
        running = min(running, p_values[idx] * n / rank)
        adjusted[idx] = running
    return adjusted


def consolidate_worker_csvs(results_dir: Path):
    """If worker subdirectories exist (e.g. workers/R004), consolidate their CSVs into tables/ and results_dir/."""
    workers_dir = results_dir / "workers"
    if not workers_dir.exists():
        return

    tables_dir = results_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    csv_names = [
        "phase_summary.csv",
        "core_members.csv",
        "core_overlap.csv",
        "ensemble_activity.csv",
        "singular_values.csv",
        "event_thresholds.csv",
        "core_shuffle_summary.csv",
        "global_factor_phase_activity.csv",
    ]

    worker_dirs = [d for d in workers_dir.iterdir() if d.is_dir()]
    if not worker_dirs:
        return

    print(f"Consolidando CSVs de {len(worker_dirs)} workers ({[w.name for w in worker_dirs]})...")
    for csv_name in csv_names:
        dfs = []
        for w in worker_dirs:
            csv_path = w / csv_name
            if csv_path.exists():
                try:
                    df = pd.read_csv(csv_path)
                    if not df.empty:
                        dfs.append(df)
                except Exception as e:
                    print(f"Aviso al leer {csv_path}: {e}")
        if dfs:
            merged_df = pd.concat(dfs, ignore_index=True)
            merged_df.drop_duplicates(inplace=True)
            merged_df.to_csv(tables_dir / csv_name, index=False)
            merged_df.to_csv(results_dir / csv_name, index=False)
            print(f"  [OK] Consolidado {csv_name}: {len(merged_df)} filas -> {tables_dir / csv_name}")


def load_data(results_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load core_overlap, phase_summary, and global_factor_phase_activity."""
    consolidate_worker_csvs(results_dir)

    overlap_file = results_dir / "tables" / "core_overlap.csv"
    if not overlap_file.exists():
        overlap_file = results_dir / "core_overlap.csv"

    summary_file = results_dir / "tables" / "phase_summary.csv"
    if not summary_file.exists():
        summary_file = results_dir / "phase_summary.csv"

    global_file = results_dir / "tables" / "global_factor_phase_activity.csv"
    if not global_file.exists():
        global_file = results_dir / "global_factor_phase_activity.csv"

    df_overlap = pd.read_csv(overlap_file) if overlap_file.exists() else pd.DataFrame()
    df_summary = pd.read_csv(summary_file) if summary_file.exists() else pd.DataFrame()
    df_global = pd.read_csv(global_file) if global_file.exists() else pd.DataFrame()

    return df_overlap, df_summary, df_global


def plot_session_jaccard_heatmap(sub_overlap: pd.DataFrame, title: str, output_path: Path):
    """Plot a clean 4x4 matrix of max Jaccard between phases for a single session."""
    is_habl = "HabL" in title
    phases = PHASES_HABL if is_habl else PHASES_TASK

    matrix = np.zeros((4, 4))
    annot = np.empty((4, 4), dtype=object)

    for i, p1 in enumerate(phases):
        matrix[i, i] = 1.0
        annot[i, i] = "1.00"
        for j, p2 in enumerate(phases):
            if i >= j:
                continue
            pair = sub_overlap[
                ((sub_overlap["phase_a"] == p1) & (sub_overlap["phase_b"] == p2)) |
                ((sub_overlap["phase_b"] == p1) & (sub_overlap["phase_a"] == p2))
            ]
            if len(pair) > 0:
                best = pair.sort_values("jaccard", ascending=False).iloc[0]
                j_val = best["jaccard"]
                p_val = best["p_permutation_unadjusted"]
                matrix[i, j] = j_val
                matrix[j, i] = j_val
                sig_star = " *" if p_val <= 0.05 else ""
                annot[i, j] = f"{j_val:.2f}{sig_star}"
                annot[j, i] = f"{j_val:.2f}{sig_star}"
            else:
                annot[i, j] = "-"
                annot[j, i] = "-"

    fig, ax = plt.subplots(figsize=(4.5, 4.0))
    im = ax.imshow(matrix, cmap="YlGnBu", vmin=0, vmax=0.8, aspect="auto")

    ax.set_xticks(range(4))
    ax.set_yticks(range(4))
    ax.set_xticklabels(phases, fontweight="bold")
    ax.set_yticklabels(phases, fontweight="bold")

    for i in range(4):
        for j in range(4):
            val = matrix[i, j]
            color = "white" if val > 0.45 else "black"
            ax.text(j, i, annot[i, j], ha="center", va="center", color=color, fontsize=10)

    ax.set_title(title, pad=12, fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Max Jaccard Core Overlap")
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def plot_global_factor_dynamics(sub_global: pd.DataFrame, title: str, output_path: Path):
    """Plot bar chart of event rate (events/min) for each global factor across the 4 phases."""
    if sub_global.empty:
        return

    factors = sorted(sub_global["ensemble"].unique())
    n_factors = len(factors)
    if n_factors == 0:
        return

    phases = PHASES_HABL if "HabL" in title else PHASES_TASK
    x = np.arange(len(phases))
    width = 0.8 / max(n_factors, 1)

    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    palette = plt.cm.tab10(np.linspace(0, 1, max(n_factors, 10)))

    for idx, f in enumerate(factors):
        f_data = sub_global[sub_global["ensemble"] == f]
        rates = []
        for p in phases:
            match = f_data[f_data["phase"] == p]
            rates.append(match["factor_vectors_per_min"].values[0] if len(match) > 0 else 0.0)

        offset = (idx - n_factors / 2 + 0.5) * width
        bars = ax.bar(x + offset, rates, width, label=f"Factor {f}", color=palette[idx], alpha=0.9)

    ax.set_xticks(x)
    ax.set_xticklabels(phases, fontweight="bold")
    ax.set_ylabel("Tasa de activacion (eventos / min)", fontweight="bold")
    ax.set_title(f"Dinamica de Factores Globales: {title}", pad=12, fontweight="bold")
    ax.legend(title="Ensamble Global", bbox_to_anchor=(1.04, 1), loc="upper left", frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def plot_consolidated_grid(df_overlap: pd.DataFrame, output_path: Path):
    """Generate the 3x3 publication master panel: HabL vs SD VEH vs SD CNO across R004, R005, R006."""
    animals = ["R004", "R005", "R006"]
    conditions = [
        ("Control (HabL)", "HabL"),
        ("SD: Control (VEH)", "VEH"),
        ("SD: Silenciamiento (CNO)", "CNO")
    ]

    fig, axes = plt.subplots(3, 3, figsize=(13, 11), sharex=False, sharey=False)

    for r, anim in enumerate(animals):
        for c, (col_title, cond_key) in enumerate(conditions):
            ax = axes[r, c]

            # filter data
            if cond_key == "HabL":
                sub = df_overlap[(df_overlap["animal"] == anim) & (df_overlap["day_name"].str.contains("HabL", case=False))]
                phases = PHASES_HABL
            else:
                sub = df_overlap[
                    (df_overlap["animal"] == anim) &
                    (df_overlap["day_name"].str.contains("SD", case=False)) &
                    (df_overlap["day_name"].str.contains(cond_key, case=False))
                ]
                phases = PHASES_TASK

            matrix = np.zeros((4, 4))
            annot = np.empty((4, 4), dtype=object)

            for i, p1 in enumerate(phases):
                matrix[i, i] = 1.0
                annot[i, i] = "1.0"
                for j, p2 in enumerate(phases):
                    if i >= j:
                        continue
                    pair = sub[
                        ((sub["phase_a"] == p1) & (sub["phase_b"] == p2)) |
                        ((sub["phase_b"] == p1) & (sub["phase_a"] == p2))
                    ]
                    if len(pair) > 0:
                        best = pair.sort_values("jaccard", ascending=False).iloc[0]
                        j_val = best["jaccard"]
                        p_val = best["p_permutation_unadjusted"]
                        matrix[i, j] = j_val
                        matrix[j, i] = j_val
                        sig_star = "*" if p_val <= 0.05 else ""
                        annot[i, j] = f"{j_val:.2f}{sig_star}"
                        annot[j, i] = f"{j_val:.2f}{sig_star}"
                    else:
                        annot[i, j] = "-"
                        annot[j, i] = "-"

            im = ax.imshow(matrix, cmap="YlGnBu", vmin=0, vmax=0.8, aspect="auto")

            ax.set_xticks(range(4))
            ax.set_yticks(range(4))
            ax.set_xticklabels(["OF1", "S", "T", "OF2"] if cond_key != "HabL" else ["1", "2", "3", "4"], fontsize=8)
            ax.set_yticklabels(["OF1", "S", "T", "OF2"] if cond_key != "HabL" else ["1", "2", "3", "4"], fontsize=8)

            for i in range(4):
                for j in range(4):
                    val = matrix[i, j]
                    color = "white" if val > 0.45 else "black"
                    ax.text(j, i, annot[i, j], ha="center", va="center", color=color, fontsize=8)

            if r == 0:
                ax.set_title(col_title, pad=10, fontweight="bold", fontsize=11)
            if c == 0:
                ax.set_ylabel(f"{anim}", fontweight="bold", fontsize=12, labelpad=10)

    fig.subplots_adjust(right=0.88, wspace=0.3, hspace=0.3)
    cbar_ax = fig.add_axes([0.91, 0.25, 0.02, 0.5])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("Similitud Jaccard de Cores (* p <= 0.05)", fontweight="bold")

    fig.suptitle("Estabilidad y Continuidad de Ensambles en CA1: HabL vs SD (VEH y CNO)", fontsize=14, fontweight="bold", y=0.98)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)
    print(f"Master panel saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate publication figures from Stoixeion results.")
    parser.add_argument("--results-dir", type=Path, default=Path("results/stoixeion/phasewise"),
                        help="Path to folder containing CSV exports.")
    parser.add_argument("--output-dir", type=Path, default=Path("results/stoixeion/publication_figures"),
                        help="Path to save generated PNG figures.")
    args = parser.parse_args()

    df_overlap, df_summary, df_global = load_data(args.results_dir)
    if df_overlap.empty:
        print(f"No core_overlap.csv found in {args.results_dir}. Aborting.")
        return

    # 1. Master Panel
    plot_consolidated_grid(df_overlap, args.output_dir / "master_panel_3x3.png")

    # 2. Individual Heatmaps
    for (anim, day), sub in df_overlap.groupby(["animal", "day_name"]):
        if ("SD" in day or "HabL" in day) and anim in ["R004", "R005", "R006"]:
            out_file = args.output_dir / "heatmaps" / f"{anim}_{day}_jaccard_heatmap.png"
            plot_session_jaccard_heatmap(sub, f"{anim} - {day}", out_file)

    # 3. Global Factor Activity Dynamics
    if not df_global.empty:
        for (anim, day), sub in df_global.groupby(["animal", "day_name"]):
            if ("SD" in day or "HabL" in day) and anim in ["R004", "R005", "R006"]:
                out_file = args.output_dir / "dynamics" / f"{anim}_{day}_global_activity.png"
                plot_global_factor_dynamics(sub, f"{anim} - {day}", out_file)

    print(f"\nTodas las figuras generadas exitosamente en: {args.output_dir}")


if __name__ == "__main__":
    main()
