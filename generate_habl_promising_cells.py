"""Generate the single canonical HabL place-cell candidate collection.

Every candidate is presented in one format:
1. trajectory with every S > 0 sample, colored by S amplitude;
2. binary S > 2 SD spatial map;
3. spatial map using every positive S amplitude.

Maps use 1-cm evaluation bins and a Gaussian sigma of 4 cm.  There are no
shuffles or split-half selection criteria in this presentation.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from core.data_loader import load_session, mapping_rows_for_columns
from core.place_fields import gaussian_rate_maps_matlab


SESSION_FILE = "2026_06_18_merged.mat"
OUTPUT_DIR = Path("results/HabL_promising_cells")
BIN_SIZE_CM = 1.0
SIGMA_CM = 4.0
THRESHOLD_SD = 2.0

# The first five are the user-selected cells. Global cell 185 is shown in OF2
# (local 128), as requested, rather than its originally selected OF4 local 142.
# The remaining cells are the visually reviewed additions; the last one was
# added in the final continuous-activity-prioritized search.
CANDIDATES = (
    ("selected", "R004", 3, 153),
    ("selected", "R004", 2, 128),
    ("selected", "R005", 1, 36),
    ("selected", "R005", 2, 108),
    ("selected", "R005", 4, 237),
    ("new", "R005", 1, 464),
    ("new", "R005", 4, 214),
    ("new", "R004", 2, 141),
    ("new", "R005", 1, 293),
    ("new", "R005", 4, 160),
    ("new", "R005", 2, 69),
    ("new", "R005", 2, 418),
    ("new", "R005", 3, 48),
    ("new", "R005", 3, 452),
    ("new", "R005", 3, 115),
    ("new", "R004", 2, 105),
    ("new", "R005", 1, 355),
    ("new", "R005", 1, 31),
)


def prepare_candidates():
    data_cache = {}
    grouped = defaultdict(list)
    records = []
    for rank, (group, animal, session, column) in enumerate(
        CANDIDATES, start=1
    ):
        data = data_cache.setdefault(animal, load_session(animal, SESSION_FILE))
        sub = data["subsessions"][session - 1]
        mapping_rows = mapping_rows_for_columns(
            data["mapping"], session - 1, sub["S"].shape[1]
        )
        record = {
            "rank": rank,
            "group": group,
            "animal": animal,
            "session": session,
            "local_cell": column,
            "global_cell": int(mapping_rows[column]),
        }
        records.append(record)
        grouped[(animal, session)].append(record)

    exact = {}
    for (animal, session), session_records in grouped.items():
        sub = data_cache[animal]["subsessions"][session - 1]
        columns = [record["local_cell"] for record in session_records]
        raw = np.asarray(sub["S"], dtype=float)[:, columns]
        x_all = np.asarray(sub["x"], dtype=float)
        y_all = np.asarray(sub["y"], dtype=float)
        t_all = np.asarray(sub["t"], dtype=float)
        valid = (
            np.isfinite(x_all)
            & np.isfinite(y_all)
            & np.isfinite(t_all)
            & np.isfinite(raw).all(axis=1)
        )
        x = x_all[valid]
        y = y_all[valid]
        t = t_all[valid]
        signals = raw[valid]
        thresholds = THRESHOLD_SD * np.std(signals, axis=0, ddof=1)
        binary = (signals > thresholds[None, :]).astype(np.float32)
        positive_signals = np.maximum(signals, 0.0)
        activity = np.column_stack((binary, positive_signals)).astype(
            np.float32
        )
        maps, edges, metadata = gaussian_rate_maps_matlab(
            x,
            y,
            t,
            activity,
            bin_size=BIN_SIZE_CM,
            sigma=SIGMA_CM,
        )
        count = len(columns)
        for index, record in enumerate(session_records):
            key = (
                record["animal"], record["session"], record["local_cell"]
            )
            exact[key] = {
                "x": x,
                "y": y,
                "signal": signals[:, index],
                "binary": binary[:, index],
                "threshold": float(thresholds[index]),
                "binary_map": maps[index],
                "full_map": maps[count + index],
                "edges": edges,
                "sample_rate_hz": metadata["sample_rate_hz"],
            }
    return records, exact


def style_axis(axis):
    axis.set_aspect("equal")
    axis.set_xlim(-25, 25)
    axis.set_ylim(-25, 25)
    axis.set_xticks((-25, 0, 25))
    axis.set_yticks((-25, 0, 25))
    axis.tick_params(labelsize=7)
    axis.set_xlabel("x (cm)", fontsize=7)
    axis.set_ylabel("y (cm)", fontsize=7)


def draw_trajectory(
    axis, record, values, include_global=False, add_colorbar=False
):
    positive = values["signal"] > 0
    positive_signal = values["signal"][positive]
    axis.plot(
        values["x"], values["y"], color="0.82", linewidth=0.32,
        rasterized=True, zorder=1,
    )
    scatter = axis.scatter(
        values["x"][positive], values["y"][positive],
        c=positive_signal, cmap="turbo", vmin=0,
        vmax=float(np.max(positive_signal)) if positive_signal.size else 1.0,
        s=5, alpha=0.82,
        linewidths=0, rasterized=True, zorder=2,
    )
    identity = (
        f" | global {record['global_cell']}" if include_global else ""
    )
    axis.set_title(
        f"#{record['rank']} {record['animal']} OF{record['session']} "
        f"local {record['local_cell']}{identity}\n"
        f"{int(positive.sum())} S > 0 samples",
        fontsize=8,
    )
    style_axis(axis)
    if add_colorbar:
        colorbar = plt.colorbar(
            scatter, ax=axis, fraction=0.046, pad=0.035
        )
        colorbar.ax.tick_params(labelsize=7)
        colorbar.set_label("S", fontsize=7)


def draw_map(axis, activity_map, edges, title, add_colorbar=False):
    peak = float(np.nanmax(activity_map))
    image = axis.imshow(
        activity_map.T,
        origin="lower",
        extent=[edges[0], edges[-1], edges[0], edges[-1]],
        cmap="turbo",
        vmin=0,
        vmax=peak if peak > 0 else 1.0,
        interpolation="nearest",
        rasterized=True,
    )
    axis.set_title(f"{title} | max {peak:.2f} act/s", fontsize=8)
    style_axis(axis)
    if add_colorbar:
        colorbar = plt.colorbar(image, ax=axis, fraction=0.046, pad=0.035)
        colorbar.ax.tick_params(labelsize=7)
        colorbar.set_label("Activity / s", fontsize=7)
    return peak


def draw_candidate_row(axes, record, values, individual=False):
    draw_trajectory(
        axes[0], record, values, include_global=not individual,
        add_colorbar=individual,
    )
    binary_peak = draw_map(
        axes[1], values["binary_map"], values["edges"],
        "S > 2 SD", add_colorbar=individual,
    )
    full_peak = draw_map(
        axes[2], values["full_map"], values["edges"],
        "All recordings S > 0", add_colorbar=individual,
    )
    return binary_peak, full_peak


def save_individuals(records, exact):
    index_rows = []
    for record in records:
        key = (record["animal"], record["session"], record["local_cell"])
        values = exact[key]
        fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.35))
        binary_peak, full_peak = draw_candidate_row(
            axes, record, values, individual=True
        )
        fig.suptitle(
            f"{record['animal']} | global cell {record['global_cell']}",
            fontsize=12, fontweight="semibold", y=0.985,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.94), w_pad=1.0)
        filename = (
            f"{record['rank']:02d}_{record['animal']}_global_"
            f"{record['global_cell']:04d}_OF{record['session']}_"
            f"local_{record['local_cell']:04d}.png"
        )
        output = OUTPUT_DIR / filename
        fig.savefig(output, dpi=220, facecolor="white")
        plt.close(fig)
        index_rows.append({
            **record,
            "positive_S_samples": int(np.sum(values["signal"] > 0)),
            "binary_spikes": int(np.sum(values["binary"])),
            "threshold_3sd": values["threshold"],
            "thresholded_max_act_per_s": binary_peak,
            "all_positive_S_max_act_per_s": full_peak,
            "bin_cm": BIN_SIZE_CM,
            "sigma_cm": SIGMA_CM,
            "png": str(output),
        })
    return index_rows


def save_contact_sheets(records, exact, per_page=6):
    outputs = []
    for page, start in enumerate(range(0, len(records), per_page), start=1):
        page_records = records[start:start + per_page]
        fig, axes = plt.subplots(
            len(page_records), 3,
            figsize=(9.2, 2.72 * len(page_records)),
            squeeze=False,
        )
        for offset, record in enumerate(page_records):
            key = (
                record["animal"], record["session"], record["local_cell"]
            )
            draw_candidate_row(axes[offset], record, exact[key])
        fig.suptitle(
            "Promising HabL place-cell candidates\n"
            "trajectory with S > 0 | S > 2 SD | all recordings S > 0",
            fontsize=13,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.965))
        output = OUTPUT_DIR / f"contact_sheet_{page:02d}.png"
        fig.savefig(output, dpi=200, facecolor="white")
        plt.close(fig)
        outputs.append(output)
    return outputs


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    records, exact = prepare_candidates()
    index_rows = save_individuals(records, exact)
    contact_sheets = save_contact_sheets(records, exact)
    index_path = OUTPUT_DIR / "candidates.csv"
    with index_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=index_rows[0].keys())
        writer.writeheader()
        writer.writerows(index_rows)
    print(
        f"Generated {len(records)} candidates and "
        f"{len(contact_sheets)} contact sheets in {OUTPUT_DIR}",
        flush=True,
    )


if __name__ == "__main__":
    main()
