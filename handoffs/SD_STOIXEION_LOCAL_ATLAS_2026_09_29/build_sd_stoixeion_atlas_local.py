"""Build a six-day, within-day Stoixeion atlas from local MAT and saved fits.

The factor fit comes from the exported MATLAB Stoixeion CSV files. This script
reconstructs the binarized population activity from the local MAT files,
validates every global assignment against that raster, recomputes phase curves
and core-cell activity, and renders readable four-phase figures. It does not
rerun the quadratic global SVD.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
from scipy.io import loadmat


PHASES = ("OF1", "SAMPLE", "TEST", "OF2")
SPEC = (
    ("R004", "VEH", "T2_SD_VEH_", "2026_06_21_merged.mat"),
    ("R004", "CNO", "T1_SD_CNO", "2026_06_19_merged.mat"),
    ("R005", "VEH", "T1_SD_VEH", "2026_06_19_merged.mat"),
    ("R005", "CNO", "T2_SD_CNO", "2026_06_21_merged.mat"),
    ("R006", "VEH", "T1_SD_VEH", "2026_06_20_merged.mat"),
    ("R006", "CNO", "T2_SD_CNO", "2026_06_22_merged.mat"),
)
PHASE_COLOR = ("#eef2f5", "#e5f2f5", "#eef2f5", "#e5f2f5")
FACTOR_COLORS = ("#16739c", "#dd6b20", "#238c73", "#8a63b0", "#b88716", "#b74662", "#596776", "#6262a8", "#7c8c2b", "#995d47")


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def mapped_phase_columns(mapping: np.ndarray, global_ids: np.ndarray, phase_index: int, n_columns: int) -> np.ndarray:
    local_ids = mapping[:, phase_index]
    rows = np.flatnonzero(np.isfinite(local_ids))
    if len(rows) != n_columns:
        raise ValueError(f"Phase {phase_index + 1}: {len(rows)} mapping IDs but {n_columns} S columns")
    ordered = rows[np.argsort(local_ids[rows])]
    if len(np.unique(local_ids[ordered])) != len(ordered):
        raise ValueError(f"Phase {phase_index + 1}: duplicate local mapping ID")
    lookup = {int(row): col for col, row in enumerate(ordered)}
    return np.array([lookup[int(row)] for row in global_ids], dtype=int)


def load_mapped_day(path: Path, global_ids: np.ndarray) -> tuple[list[np.ndarray], list[np.ndarray], np.ndarray]:
    act = loadmat(path, simplify_cells=True, variable_names=["act"])["act"]
    mapping = np.asarray(act["mapping"], dtype=float)
    if mapping.shape[1] != 4:
        raise ValueError(f"Expected 4 mapping columns in {path}")
    common = set(np.flatnonzero(np.isfinite(mapping).all(axis=1)).tolist())
    if set(global_ids.tolist()) != common:
        raise ValueError(f"Exported global IDs do not equal 4/4 MAT mapping IDs in {path}")
    signals: list[np.ndarray] = []
    times: list[np.ndarray] = []
    for i in range(4):
        source = np.asarray(act["S"][i])
        columns = mapped_phase_columns(mapping, global_ids, i, source.shape[1])
        signals.append(np.asarray(source[:, columns], dtype=np.float32))
        t = np.asarray(act["t"][i], dtype=float).ravel()
        if len(t) != len(source):
            raise ValueError(f"Time/S row mismatch for {path}, phase {i + 1}")
        times.append(t - t[0])
    del act
    return signals, times, mapping


def duration_s(t: np.ndarray) -> float:
    dt = np.diff(t)
    dt = dt[np.isfinite(dt) & (dt > 0)]
    return float(t[-1] + (np.median(dt) if len(dt) else 0.05))


def binned_curve(t: np.ndarray, selected: np.ndarray, significant: np.ndarray, start_s: float, phase: str, factor: int) -> list[dict]:
    stop = duration_s(t)
    edges = np.arange(0.0, np.ceil(stop / 10.0) * 10.0 + 10.01, 10.0)
    all_counts = np.histogram(t, bins=edges)[0]
    selected_counts = np.histogram(t[selected], bins=edges)[0]
    sig_counts = np.histogram(t[significant], bins=edges)[0]
    records = []
    for j, total in enumerate(all_counts):
        if total == 0:
            continue
        records.append({
            "phase": phase, "factor": factor, "bin_start_s": float(edges[j]),
            "bin_center_day_s": float(start_s + (edges[j] + edges[j + 1]) / 2.0),
            "n_frames": int(total), "n_significant": int(sig_counts[j]),
            "n_assigned": int(selected_counts[j]),
            "pct_all_frames": float(100.0 * selected_counts[j] / total),
            "pct_significant_frames": float(100.0 * selected_counts[j] / sig_counts[j]) if sig_counts[j] else float("nan"),
        })
    return records


def count_bouts(local_indices: np.ndarray) -> int:
    if not len(local_indices):
        return 0
    return int(1 + np.count_nonzero(np.diff(np.sort(local_indices)) > 1))


def binary_cosine_sample(events: np.ndarray, significant_global: np.ndarray, phase_bounds: np.ndarray, scut: float, per_phase: int = 110) -> tuple[np.ndarray, np.ndarray, list[int]]:
    # TF-IDF uses the full set of selected high-activity vectors. Only the
    # displayed pairwise matrix is sampled, equally across the four phases.
    high = events[significant_global]
    df = high.sum(axis=0).astype(float)
    idf = 1.0 + np.log(len(significant_global) / np.maximum(df, 1.0))
    chosen = []
    counts = []
    for i in range(4):
        frames = significant_global[(significant_global >= phase_bounds[i]) & (significant_global < phase_bounds[i + 1])]
        n = min(per_phase, len(frames))
        sample = frames[np.linspace(0, len(frames) - 1, n, dtype=int)] if n else np.array([], dtype=int)
        chosen.append(sample)
        counts.append(len(sample))
    selected = np.concatenate(chosen)
    x = events[selected].astype(np.float32) * idf.astype(np.float32)[None, :]
    x /= np.maximum(np.linalg.norm(x, axis=1, keepdims=True), 1e-12)
    sim = x @ x.T
    return (sim > scut), selected, counts


def plot_timeline(path: Path, title: str, factors: list[int], curves: list[dict], assignments: list[dict], bounds_s: np.ndarray, catalog: list[dict]) -> None:
    n = len(factors)
    fig, axes = plt.subplots(n, 1, figsize=(17.5, max(7, 2.55 * n + 1.6)), sharex=True, layout="constrained", squeeze=False)
    axes = axes[:, 0]
    cat = {int(r["factor"]): r for r in catalog}
    for ax, factor in zip(axes, factors):
        color = FACTOR_COLORS[(factor - 1) % len(FACTOR_COLORS)]
        sub = sorted((r for r in curves if r["factor"] == factor), key=lambda r: r["bin_center_day_s"])
        peak = max([r["pct_all_frames"] for r in sub], default=0.0)
        ymax = max(1.0, peak * 1.22)
        for i in range(4):
            ax.axvspan(bounds_s[i] / 60, bounds_s[i + 1] / 60, color=PHASE_COLOR[i], lw=0, zorder=0)
        for i, phase in enumerate(PHASES):
            segment = [r for r in sub if r["phase"] == phase]
            if segment:
                x = np.array([r["bin_center_day_s"] / 60 for r in segment])
                y = np.array([r["pct_all_frames"] for r in segment])
                ax.fill_between(x, y, color=color, alpha=0.23, step="mid", lw=0)
                ax.plot(x, y, color=color, lw=1.8, drawstyle="steps-mid")
            if i:
                ax.axvline(bounds_s[i] / 60, color="#a9b9c0", lw=0.9)
        ticks = np.array([r["day_elapsed_s"] for r in assignments if r["factor"] == factor], dtype=float) / 60
        if len(ticks):
            ax.scatter(ticks, np.full(len(ticks), ymax * 0.97), marker="|", s=17, linewidths=0.55, color=color, alpha=0.55, rasterized=True)
        ax.set_ylim(0, ymax * 1.06)
        ax.set_yticks([0, round(ymax / 2, 1), round(ymax, 1)])
        ax.tick_params(axis="y", labelsize=8)
        ax.set_ylabel(f"Factor {factor}\n% frames", rotation=0, ha="right", va="center", labelpad=30, fontweight="bold", fontsize=10)
        ax.text(0.995, 0.78, f"{int(cat[factor]['n_core'])} core · {int(cat[factor]['n_members'])} miembros", transform=ax.transAxes, ha="right", va="center", fontsize=9, color="#455a64", bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"})
        ax.grid(axis="y", color="#dce5e9", lw=0.6)
        ax.spines[["top", "right"]].set_visible(False)
    for i, phase in enumerate(PHASES):
        axes[0].text((bounds_s[i] + bounds_s[i + 1]) / (2 * bounds_s[-1]), 1.08, phase, transform=axes[0].transAxes, ha="center", va="bottom", fontweight="bold", color="#1e5067", fontsize=11)
    axes[-1].set_xlim(0, bounds_s[-1] / 60)
    axes[-1].set_xlabel("Minutos desde el inicio de OF1", fontsize=11)
    fig.suptitle(title + " · actividad de todos los factores globales", x=0.06, ha="left", fontsize=17, fontweight="bold")
    fig.savefig(path, dpi=175, facecolor="white")
    plt.close(fig)


def plot_method(path: Path, title: str, events: np.ndarray, significant: np.ndarray, binary_m: np.ndarray, m_counts: list[int], sv: pd.DataFrame, factors: list[int], curves: list[dict], phase_bounds: np.ndarray, bounds_s: np.ndarray, scut: float) -> None:
    fig = plt.figure(figsize=(17.5, 12), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[1.05, 0.95], width_ratios=[1.1, 1.0])
    ax_raster = fig.add_subplot(grid[0, 0])
    ax_m = fig.add_subplot(grid[0, 1])
    ax_svd = fig.add_subplot(grid[1, 0])
    ax_heat = fig.add_subplot(grid[1, 1])

    # Pool roughly 1 s of frames per image column: preserves an interpretable
    # whole-day raster while making the downsampling explicit.
    raster_parts = []
    raster_lengths = []
    for i in range(4):
        sub = events[phase_bounds[i]:phase_bounds[i + 1]]
        pooled = np.stack([sub[j:j + 20].any(axis=0) for j in range(0, len(sub), 20)], axis=0)
        raster_parts.append(pooled)
        raster_lengths.append(len(pooled))
    raster = np.concatenate(raster_parts).T
    ax_raster.imshow(raster, cmap=ListedColormap(["white", "#263e4b"]), vmin=0, vmax=1, aspect="auto", interpolation="nearest", rasterized=True)
    for x in np.cumsum(raster_lengths)[:-1]:
        ax_raster.axvline(x - 0.5, color="#ef9a5b", lw=1)
    ax_raster.set_title("A · actividad binaria mapeada (máximo por ≈1 s)", loc="left", fontsize=12, fontweight="bold")
    ax_raster.set_ylabel("Célula presente en las 4 fases")
    ax_raster.set_xticks((np.cumsum([0] + raster_lengths)[:-1] + np.cumsum(raster_lengths)) / 2, PHASES)

    ax_m.imshow(binary_m, cmap=ListedColormap(["white", "#263e4b"]), vmin=0, vmax=1, aspect="equal", interpolation="nearest", rasterized=True)
    edges = np.cumsum([0] + m_counts)
    for x in edges[1:-1]:
        ax_m.axvline(x - 0.5, color="#ef9a5b", lw=0.8)
        ax_m.axhline(x - 0.5, color="#ef9a5b", lw=0.8)
    centers = (edges[:-1] + edges[1:]) / 2
    ax_m.set_xticks(centers, PHASES)
    ax_m.set_yticks(centers, PHASES)
    ax_m.set_title(f"B · coseno TF-IDF > scut {scut:.2f} (muestra)", loc="left", fontsize=12, fontweight="bold")

    sv = sv.sort_values("singular_rank")
    selected = sv[sv.selected_by_stoixeion == 1]
    shown = max(30, min(80, 8 * len(selected)))
    sv_shown = sv[sv.singular_rank <= shown]
    ax_svd.plot(sv_shown.singular_rank, sv_shown.singular_value, color="#263e4b", lw=1.6)
    ax_svd.scatter(selected.singular_rank, selected.singular_value, color="#dd6b20", s=34, label="Seleccionado", zorder=3)
    ax_svd.set_yscale("log")
    ax_svd.set_xlim(0, shown)
    ax_svd.set_xlabel("Rango singular")
    ax_svd.set_ylabel("Valor singular")
    ax_svd.set_title(f"C · primeros {shown} valores singulares exportados", loc="left", fontsize=12, fontweight="bold")
    ax_svd.legend(loc="upper right", frameon=False, fontsize=9)
    ax_svd.spines[["top", "right"]].set_visible(False)

    bins = sorted({float(r["bin_center_day_s"]) for r in curves})
    b_index = {v: i for i, v in enumerate(bins)}
    values = np.zeros((len(factors), len(bins)), dtype=float)
    for r in curves:
        values[factors.index(int(r["factor"])), b_index[float(r["bin_center_day_s"])]] = float(r["pct_all_frames"])
    vmax = max(1.0, float(np.quantile(values, 0.99)))
    im = ax_heat.imshow(values, cmap="magma", vmin=0, vmax=vmax, aspect="auto", interpolation="nearest")
    for x in bounds_s[1:-1]:
        idx = int(np.searchsorted(bins, x))
        ax_heat.axvline(idx - 0.5, color="white", lw=0.9, alpha=0.9)
    centers = [int(np.searchsorted(bins, (bounds_s[i] + bounds_s[i + 1]) / 2)) for i in range(4)]
    ax_heat.set_xticks(centers, PHASES)
    ax_heat.set_yticks(range(len(factors)), [f"Factor {f}" for f in factors])
    ax_heat.set_title("D · frames asignados a cada factor (% por 10 s)", loc="left", fontsize=12, fontweight="bold")
    fig.colorbar(im, ax=ax_heat, label="% de todos los frames; color saturado al p99", shrink=0.75)
    fig.suptitle(title + " · lectura tipo Carrillo-Reid", x=0.02, ha="left", fontsize=17, fontweight="bold")
    fig.savefig(path, dpi=170, facecolor="white")
    plt.close(fig)


def process_day(spec: tuple[str, str, str, str], data_dir: Path, exports: Path, out: Path, tables: dict[str, pd.DataFrame]) -> dict:
    animal, treatment, day_name, session = spec
    key = (tables["summary"].animal == animal) & (tables["summary"].day_name == day_name) & (tables["summary"].session == session)
    summary = tables["summary"][key & (tables["summary"].phase == "GLOBAL_CONCAT")]
    if len(summary) != 1 or summary.iloc[0].status != "ok":
        raise ValueError(f"No successful global fit for {animal} {day_name}")
    fit = summary.iloc[0]
    dkey = lambda df: (df.animal == animal) & (df.day_name == day_name) & (df.session == session)
    threshold_rows = tables["thresholds"][dkey(tables["thresholds"]) & (tables["thresholds"].phase == "GLOBAL_CONCAT")].sort_values("local_cell_index")
    if len(threshold_rows) != int(fit.n_common_cells):
        raise ValueError(f"Threshold count mismatch for {animal} {day_name}")
    global_ids = threshold_rows.mapped_cell_id.to_numpy(dtype=int)
    if len(np.unique(global_ids)) != len(global_ids):
        raise ValueError("Duplicate mapped cell ID in exported thresholds")
    thresholds = threshold_rows.threshold.to_numpy(dtype=np.float32)
    mat_path = data_dir / animal / session
    signals, times, mapping = load_mapped_day(mat_path, global_ids)
    phase_n = [len(s) for s in signals]
    phase_bounds = np.cumsum([0] + phase_n)
    phase_durations = [duration_s(t) for t in times]
    bounds_s = np.cumsum([0.0] + phase_durations)
    S = np.concatenate(signals)
    recalculated = 3 * np.std(S.astype(np.float64), axis=0, ddof=1)
    threshold_difference = float(np.max(np.abs(recalculated - thresholds)))
    # MATLAB computed the thresholds on the same concatenated mapped traces.
    if threshold_difference > 1e-4:
        raise ValueError(f"Exported threshold mismatch ({threshold_difference:g}) for {animal} {day_name}")
    events = np.isfinite(S) & (S > thresholds[None, :])
    pks = int(fit.pks)
    significant = np.flatnonzero(events.sum(axis=1) >= pks)
    if len(significant) != int(fit.n_significant_vectors):
        raise ValueError(f"High-activity vector mismatch {animal} {day_name}: {len(significant)} vs {int(fit.n_significant_vectors)}")
    gpa = tables["gpa"][dkey(tables["gpa"])].copy()
    for i, phase in enumerate(PHASES):
        rows = gpa[gpa.phase == phase]
        if len(rows) != int(fit.n_ensembles):
            raise ValueError(f"Factor-phase summary incomplete for {animal} {day_name} {phase}")
        n_sig = np.count_nonzero((significant >= phase_bounds[i]) & (significant < phase_bounds[i + 1]))
        if not (rows.n_significant_vectors.to_numpy(dtype=int) == n_sig).all():
            raise ValueError(f"Phase significant count mismatch for {animal} {day_name} {phase}")
        if not (rows.n_phase_frames.to_numpy(dtype=int) == phase_n[i]).all():
            raise ValueError(f"Phase frame count mismatch for {animal} {day_name} {phase}")

    member_rows = tables["members"][dkey(tables["members"]) & (tables["members"].phase == "GLOBAL_CONCAT")]
    core_ids_by_factor = {}
    catalog = []
    factors = list(range(1, int(fit.n_ensembles) + 1))
    for factor in factors:
        members = member_rows[member_rows.ensemble == factor]
        all_ids = sorted(members.mapped_cell_id.astype(int).unique().tolist())
        core_ids = sorted(members[members.is_core == 1].mapped_cell_id.astype(int).unique().tolist())
        if not all_ids or not set(core_ids).issubset(set(all_ids)) or not set(all_ids).issubset(set(global_ids.tolist())):
            raise ValueError(f"Member/core mapping mismatch {animal} {day_name} factor {factor}")
        core_ids_by_factor[factor] = core_ids
        catalog.append({"animal": animal, "treatment": treatment, "day_name": day_name, "session": session,
                        "factor": factor, "n_members": len(all_ids), "n_core": len(core_ids),
                        "member_ids_0based": ";".join(map(str, all_ids)), "core_ids_0based": ";".join(map(str, core_ids))})
    if not (set(member_rows.ensemble.astype(int)) == set(factors)):
        raise ValueError(f"Missing factor membership in {animal} {day_name}")

    raw_assignments = tables["activity"][dkey(tables["activity"]) & tables["activity"].phase.str.startswith("GLOBAL_")].copy()
    if raw_assignments.frame.duplicated().any():
        raise ValueError(f"Global frame assigned to more than one factor in {animal} {day_name}")
    assignment_rows = []
    factor_local: dict[tuple[int, str], np.ndarray] = {}
    for i, phase in enumerate(PHASES):
        sub = raw_assignments[raw_assignments.phase == "GLOBAL_" + phase]
        global_index = sub.frame.to_numpy(dtype=int) - 1
        local_index = global_index - phase_bounds[i]
        if len(local_index) and (local_index.min() < 0 or local_index.max() >= phase_n[i]):
            raise ValueError(f"Assignment outside {phase} in {animal} {day_name}")
        if len(local_index) and not np.all(events.sum(axis=1)[global_index] >= pks):
            raise ValueError(f"Assignment on non-significant frame in {animal} {day_name} {phase}")
        if len(local_index) and np.max(np.abs(times[i][local_index] - sub.time_s.to_numpy(dtype=float))) > 0.055:
            raise ValueError(f"Exported local time mismatch for {animal} {day_name} {phase}")
        for j, (_, row) in enumerate(sub.iterrows()):
            assignment_rows.append({"animal": animal, "treatment": treatment, "day_name": day_name, "session": session,
                                    "phase": phase, "factor": int(row.ensemble), "global_frame_1based": int(row.frame),
                                    "local_frame_1based": int(local_index[j] + 1), "time_s": float(times[i][local_index[j]]),
                                    "day_elapsed_s": float(bounds_s[i] + times[i][local_index[j]])})
        for factor in factors:
            factor_local[(factor, phase)] = local_index[sub.ensemble.to_numpy(dtype=int) == factor]
            reported = gpa[(gpa.phase == phase) & (gpa.ensemble == factor)].iloc[0]
            if len(factor_local[(factor, phase)]) != int(reported.factor_vectors):
                raise ValueError(f"Factor count mismatch in {animal} {day_name} {phase} F{factor}")
    if len(assignment_rows) != sum(int(row["factor_vectors"]) for _, row in gpa.iterrows()):
        raise ValueError(f"Total factor count mismatch in {animal} {day_name}")

    curves = []
    phase_stats = []
    id_to_col = {int(cid): j for j, cid in enumerate(global_ids)}
    for i, phase in enumerate(PHASES):
        sig_local = significant[(significant >= phase_bounds[i]) & (significant < phase_bounds[i + 1])] - phase_bounds[i]
        for factor in factors:
            selected = factor_local[(factor, phase)]
            curves.extend(binned_curve(times[i], selected, sig_local, bounds_s[i], phase, factor))
            cores = core_ids_by_factor[factor]
            if cores:
                corecols = [id_to_col[c] for c in cores]
                mean_core_fraction = float(events[phase_bounds[i]:phase_bounds[i + 1], corecols].mean())
            else:
                mean_core_fraction = float("nan")
            phase_stats.append({"animal": animal, "treatment": treatment, "day_name": day_name, "phase": phase,
                                "factor": factor, "n_phase_frames": phase_n[i], "n_significant_frames": len(sig_local),
                                "n_assigned_frames": len(selected), "n_assigned_bouts": count_bouts(selected),
                                "assigned_frames_per_min": float(len(selected) * 60 / phase_durations[i]),
                                "pct_all_frames": float(100 * len(selected) / phase_n[i]),
                                "pct_significant_frames": float(100 * len(selected) / len(sig_local)) if len(sig_local) else float("nan"),
                                "mean_core_cell_active_fraction": mean_core_fraction, "n_core": len(cores)})

    binary_m, chosen, m_counts = binary_cosine_sample(events, significant, phase_bounds, float(fit.scut))
    sv = tables["sv"][dkey(tables["sv"]) & (tables["sv"].phase == "GLOBAL_CONCAT")]
    if sv.empty:
        raise ValueError(f"Missing global singular values for {animal} {day_name}")
    stem = f"{animal}_{treatment}_{session[:10]}"
    title = f"{animal} · S–D · {treatment} · {session[:10].replace('_', '/')}"
    day_dir = out / "days" / stem
    day_dir.mkdir(parents=True, exist_ok=True)
    write_csv(day_dir / "factors.csv", catalog)
    write_csv(day_dir / "phase_metrics.csv", phase_stats)
    quality_rows = []
    for item in catalog:
        factor = int(item["factor"])
        by_phase = {row["phase"]: row for row in phase_stats if row["factor"] == factor}
        assigned_total = sum(by_phase[p]["n_assigned_frames"] for p in PHASES)
        dominant = max(PHASES, key=lambda p: by_phase[p]["n_assigned_frames"])
        flags = []
        if item["n_core"] < 3:
            flags.append("core_menor_3")
        if item["n_members"] / len(global_ids) > 0.8:
            flags.append("miembros_mas_80pct_celulas")
        if assigned_total / len(events) < 0.01:
            flags.append("asignaciones_menos_1pct_frames")
        quality_rows.append({"animal": animal, "treatment": treatment, "day_name": day_name,
                             "factor": factor, "n_core": item["n_core"], "n_members": item["n_members"],
                             "member_fraction_of_mapped_cells": item["n_members"] / len(global_ids),
                             "assigned_frames_total": assigned_total,
                             "assigned_fraction_of_recording": assigned_total / len(events),
                             "dominant_phase_by_count": dominant,
                             "dominant_phase_share_of_assignments": by_phase[dominant]["n_assigned_frames"] / assigned_total if assigned_total else float("nan"),
                             "n_phases_with_assignments": sum(by_phase[p]["n_assigned_frames"] > 0 for p in PHASES),
                             "descriptive_review_flags": ";".join(flags)})
    write_csv(day_dir / "factor_quality.csv", quality_rows)
    write_csv(day_dir / "frame_assignments.csv", assignment_rows)
    write_csv(day_dir / "curves_10s.csv", [{"animal": animal, "treatment": treatment, "day_name": day_name, **r} for r in curves])
    np.savez_compressed(day_dir / "cosine_binary_sample.npz", binary_similarity=binary_m.astype(np.uint8),
                        selected_global_frame_1based=chosen + 1, per_phase_sample_count=np.array(m_counts),
                        pks=np.array(pks), scut=np.array(float(fit.scut)))
    plot_timeline(day_dir / "all_ensembles_four_phases.png", title, factors, curves, assignment_rows, bounds_s, catalog)
    plot_method(day_dir / "carrillo_method_four_phases.png", title, events, significant, binary_m, m_counts, sv,
                factors, curves, phase_bounds, bounds_s, float(fit.scut))
    qc = {"animal": animal, "treatment": treatment, "day_name": day_name, "session": session,
          "n_cells_4of4": len(global_ids), "phase_frames": dict(zip(PHASES, phase_n)),
          "n_global_factors": len(factors), "n_significant_vectors_reconstructed": len(significant),
          "n_significant_vectors_exported": int(fit.n_significant_vectors),
          "n_assigned_frames": len(assignment_rows), "max_abs_threshold_difference": threshold_difference,
          "pks": pks, "scut": float(fit.scut), "mat_sha256": sha256(mat_path),
          "relative_day_dir": str(day_dir.relative_to(out)).replace("\\", "/")}
    (day_dir / "QC.json").write_text(json.dumps(qc, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"{animal} {treatment}: {len(global_ids)} cells, {len(factors)} factors, {len(assignment_rows)} assigned frames; QC OK", flush=True)
    return qc


def write_index(out: Path, qcs: list[dict]) -> None:
    cards = []
    for row in qcs:
        rel = row["relative_day_dir"]
        name = f"{row['animal']} · {row['treatment']} · {row['session'][:10]}"
        cards.append(f'''<section class="card"><h2>{html.escape(name)}</h2><p>{row['n_cells_4of4']} células 4/4 · {row['n_global_factors']} factores · {row['n_assigned_frames']:,} frames asignados</p>
        <a href="{rel}/all_ensembles_four_phases.png"><img src="{rel}/all_ensembles_four_phases.png" alt="Actividad de todos los factores"></a>
        <nav><a href="{rel}/all_ensembles_four_phases.png">Líneas de tiempo</a><a href="{rel}/carrillo_method_four_phases.png">Panel Carrillo</a><a href="{rel}/factors.csv">Células miembro y core</a><a href="{rel}/factor_quality.csv">QC descriptivo</a><a href="{rel}/phase_metrics.csv">Métricas por fase</a><a href="{rel}/frame_assignments.csv">Frames asignados</a><a href="{rel}/cosine_binary_sample.npz">Coseno binario</a></nav></section>''')
    text = '''<!doctype html><html lang="es"><meta charset="utf-8"><title>Atlas Stoixeion S–D local</title><style>
    body{font:16px/1.5 system-ui,sans-serif;color:#172c38;background:#f4f7f8;margin:0 auto;max-width:1320px;padding:2rem}
    h1{font-size:2rem;margin-bottom:.3rem}header p{max-width:950px;color:#425a65}nav{display:flex;flex-wrap:wrap;gap:.65rem;margin-top:.7rem}
    a{color:#075c83}.card{background:white;border:1px solid #d7e1e5;border-radius:14px;padding:1.25rem;margin:1.6rem 0;box-shadow:0 2px 8px #17313e10}
    .card img{width:100%;height:360px;object-fit:cover;object-position:top left;border:1px solid #e3e9ec}
    .card nav a{padding:.35rem .65rem;background:#eaf4f7;border-radius:7px;text-decoration:none}
    </style><header><h1>Stoixeion · seis jornadas S–D</h1><p>Tres animales, un día VEH y uno CNO por animal. Cada factor se ajustó una vez sobre las cuatro fases del día. Las líneas de tiempo usan los frames asignados por el ajuste MATLAB exportado, verificados aquí contra los MAT y mappings locales. Las curvas y figuras se calcularon localmente. Los factores <strong>no</strong> tienen identidad entre fechas.</p><p><a href="comparacion_seis_jornadas.png">Comparación de los seis días</a> · <a href="RESUMEN_CORTO.md">Resumen</a> · <a href="LEEME.md">Método y límites</a> · <a href="QC_ALL_DAYS.json">QC</a></p></header>''' + "\n".join(cards) + "</html>"
    (out / "index.html").write_text(text, encoding="utf-8")


def write_short_summary(out: Path, qcs: list[dict]) -> None:
    quality = pd.concat([pd.read_csv(out / q["relative_day_dir"] / "factor_quality.csv") for q in qcs], ignore_index=True)
    metrics = pd.concat([pd.read_csv(out / q["relative_day_dir"] / "phase_metrics.csv") for q in qcs], ignore_index=True)
    lines = [
        "# Seis días S–D: resultado del atlas local", "",
        f"Se reconstruyeron **{len(quality)} factores globales** de **seis jornadas** (tres animales × VEH/CNO). Los seis `.mat` pasaron el QC de mapping, umbral, número de vectores significativos y frames asignados. El ajuste SVD de los factores proviene del export MATLAB existente; las curvas, matrices muestreadas y tablas se calcularon aquí desde los `.mat`.", "",
        "| Animal | Tratamiento | Células 4/4 | Factores | Core <3 células | Frames asignados |", "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for q in qcs:
        rows = quality[(quality.animal == q["animal"]) & (quality.treatment == q["treatment"])]
        lines.append(f"| {q['animal']} | {q['treatment']} | {q['n_cells_4of4']} | {q['n_global_factors']} | {int((rows.n_core < 3).sum())} | {q['n_assigned_frames']} |")
    broad = int((quality.member_fraction_of_mapped_cells > 0.8).sum())
    small = int((quality.n_core < 3).sum())
    lines.extend(["", f"**Control de interpretabilidad:** {broad}/{len(quality)} factores incluyen como «miembros» a más del 80 % de las células mapeadas; {small}/{len(quality)} tienen un core de una o dos células. Esas banderas son descriptivas, no un filtro estadístico. La lista *core* resulta más informativa que la membresía amplia para comparar grupos celulares. En R004 CNO, los seis factores tienen core de una o dos células, por lo que no los mostraría como ensambles de muchas neuronas sincronizadas.", ""])
    r5 = metrics[(metrics.animal == "R005") & (metrics.treatment == "VEH") & (metrics.factor == 5)].set_index("phase")
    counts = [int(r5.loc[p, "n_assigned_frames"]) for p in PHASES]
    lines.extend([f"**Ejemplo destacado, R005 VEH Factor 5:** {counts[0]}/{counts[1]}/{counts[2]}/{counts[3]} frames asignados en OF1/SAMPLE/TEST/OF2. Es un patrón de ese día concentrado en TEST. El día CNO también presenta actividad tardía en TEST, pero sus factores se ajustaron de nuevo y no comparten identidad garantizada con VEH.", "", "**Conclusión:** el atlas resuelve qué factores se detectaron y cuándo se expresan dentro de cada jornada. No identifica qué información espacial representan ni demuestra un efecto CNO. Para la metodología, mirar el panel de coseno, espectro y asignaciones; para la dinámica, las figuras altas de cada día.", ""])
    (out / "RESUMEN_CORTO.md").write_text("\n".join(lines), encoding="utf-8")


def plot_six_day_overview(out: Path, qcs: list[dict]) -> None:
    lookup = {(q["animal"], q["treatment"]): q for q in qcs}
    fig, axes = plt.subplots(3, 2, figsize=(19, 13), layout="constrained")
    for row_index, animal in enumerate(("R004", "R005", "R006")):
        panels = []
        samples = []
        for treatment in ("VEH", "CNO"):
            q = lookup[(animal, treatment)]
            base = out / q["relative_day_dir"]
            curve = pd.read_csv(base / "curves_10s.csv")
            quality = pd.read_csv(base / "factor_quality.csv").set_index("factor")
            n_factors = int(q["n_global_factors"])
            data_parts = []
            lengths = []
            for phase in PHASES:
                part = curve[curve.phase == phase]
                bins = sorted(part.bin_start_s.unique())
                block = np.zeros((n_factors, len(bins)), dtype=float)
                for f in range(1, n_factors + 1):
                    vals = part[part.factor == f].sort_values("bin_start_s").pct_all_frames.to_numpy(dtype=float)
                    if len(vals) != len(bins):
                        raise ValueError(f"Incomplete overview bins: {animal} {treatment} {phase} factor {f}")
                    block[f - 1] = vals
                data_parts.append(block)
                lengths.append(len(bins))
            panels.append((treatment, np.concatenate(data_parts, axis=1), lengths, quality, q))
            samples.append(panels[-1][1])
        all_values = np.concatenate([x.ravel() for x in samples])
        vmax = max(1.0, float(np.quantile(all_values, 0.995)))
        for col_index, (treatment, matrix, lengths, quality, q) in enumerate(panels):
            ax = axes[row_index, col_index]
            im = ax.imshow(matrix, cmap="magma", vmin=0, vmax=vmax, aspect="auto", interpolation="nearest")
            borders = np.cumsum([0] + lengths)
            centers = (borders[:-1] + borders[1:]) / 2
            ax.set_xticks(centers, PHASES, fontsize=10)
            ax.set_yticks(range(len(matrix)), [f"F{f} · core {int(quality.loc[f, 'n_core'])}" for f in range(1, len(matrix) + 1)], fontsize=9)
            for boundary in borders[1:-1]:
                ax.axvline(boundary - 0.5, color="white", lw=1, alpha=0.8)
            ax.set_title(f"{animal} · {treatment} · {q['session'][:10].replace('_', '/')} · {len(matrix)} factores", fontsize=13, fontweight="bold")
        fig.colorbar(im, ax=list(axes[row_index]), shrink=0.82, pad=0.01, label=f"{animal}: % de frames / 10 s (color hasta p99,5)")
    fig.suptitle("Actividad de todos los factores Stoixeion · 3 animales × VEH/CNO · tarea S–D", fontsize=19, fontweight="bold")
    fig.savefig(out / "comparacion_seis_jornadas.png", dpi=175, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--exports-dir", type=Path, default=Path("results/stoixeion_fast_clean"))
    parser.add_argument("--out", type=Path, default=Path("results/SD_stoixeion_atlas_local_20260929"))
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    exports = args.exports_dir.resolve()
    tables = {
        "summary": pd.read_csv(exports / "phase_summary.csv"),
        "thresholds": pd.read_csv(exports / "event_thresholds.csv"),
        "members": pd.concat([pd.read_csv(exports / "workers" / animal / "ensemble_members.csv") for animal in ("R004", "R005", "R006")], ignore_index=True),
        "activity": pd.read_csv(exports / "ensemble_activity.csv"),
        "gpa": pd.read_csv(exports / "global_factor_phase_activity.csv"),
        "sv": pd.read_csv(exports / "singular_values.csv"),
    }
    qcs = [process_day(spec, args.data_dir.resolve(), exports, out, tables) for spec in SPEC]
    (out / "QC_ALL_DAYS.json").write_text(json.dumps(qcs, indent=2, ensure_ascii=False), encoding="utf-8")
    write_index(out, qcs)
    write_short_summary(out, qcs)
    plot_six_day_overview(out, qcs)
    readme = """# Atlas Stoixeion S–D reconstruido localmente

Abrí `index.html`. La lámina `comparacion_seis_jornadas.png` pone VEH y CNO lado a lado para cada animal, con una escala de color compartida **dentro de ese animal**. Además hay una figura alta con **todos** los factores globales de cada día, un panel tipo Carrillo-Reid, las listas de células miembro y core, los frames asignados y la matriz binaria de coseno muestreada. `factor_quality.csv` añade banderas **descriptivas** para cores de menos de tres células, membresía mayor al 80 % de las células mapeadas y asignaciones en menos del 1 % de los frames; no descarta factores ni constituye una prueba estadística.

## Procedencia y cálculo

- **Identidad de los factores y frames asignados:** ajuste MATLAB Stoixeion ya exportado en `results/stoixeion_fast_clean/`; cuatro fases concatenadas por jornada. Este script **no volvió a ajustar la SVD global**: para R004 CNO la matriz de 32.877 vectores ocuparía más de 8 GB solo como arreglo `float64`, antes de las copias y SVD, mientras que esta PC tiene 8 GB de RAM y no tiene MATLAB.
- **Cálculo local nuevo:** lectura de los seis `.mat` de `data/`, mapping de las células presentes en las cuatro fases, umbrales de `S` exportados (`3×SD(S)` global por célula), reconstrucción de los vectores de alta actividad, validación de cada frame asignado y de los conteos por factor/fase, curvas en ventanas de 10 s, métricas de células core y matriz de coseno TF-IDF binarizada de una muestra temporal equilibrada entre fases. `QC_ALL_DAYS.json` comprueba las seis jornadas y registra SHA-256 de cada `.mat`.
- **Matriz M:** hasta 110 vectores de alta actividad por fase; TF-IDF usa todos los vectores significativos del día, `scut` el umbral del ajuste exportado. La matriz se muestra **antes** de los dos filtros Jaccard y la SVD. No es la matriz completa.
- **Eje temporal:** `OF1 → SAMPLE → TEST → OF2`; cada raya es un frame asignado; la curva muestra porcentaje de **todos** los frames asignados en cada ventana de 10 s. La escala vertical cambia por factor y está impresa en cada fila. `phase_metrics.csv` permite comparar valores exactos entre fases.
- **IDs de células:** `mapped_cell_id` del `.mat`, base cero, válido entre fases del mismo día. No identifica la misma célula entre días VEH/CNO.

## Lectura científica

Esto permite ver cuándo se expresa cada factor a lo largo de una jornada y descargar la lista de sus células. Un factor global fue definido usando **las cuatro fases**; su presencia en SAMPLE y TEST no es una prueba independiente de reactivación. Cada marca es un frame, no un episodio independiente. Los paneles describen seis jornadas anidadas en tres animales y no demuestran por sí solos un efecto de CNO, codificación del objeto desplazado ni una explicación del comportamiento.

## Reproducción

Desde la raíz de `revision-datos`:

```powershell
python coactivation_analysis/build_sd_stoixeion_atlas_local.py --data-dir data --exports-dir results/stoixeion_fast_clean --out results/SD_stoixeion_atlas_local_20260929
```
"""
    (out / "LEEME.md").write_text(readme, encoding="utf-8")
    print(f"Atlas ready: {out / 'index.html'}", flush=True)


if __name__ == "__main__":
    main()
