"""Generate multi-session evolution figures for all cells present in >=2 HabL sessions.

Method used:
- Trajectory with every S > 0 sample colored by S amplitude.
- "All recordings" rate map using the continuous Gaussian kernel estimator
  (core.place_fields.gaussian_rate_maps_matlab) with bin_size = 1.0 cm,
  sigma = 4.0 cm, and positive continuous deconvolution amplitudes (S > 0).
- Figures display the 4 open field sessions (OF1, OF2, OF3, OF4) aligned in
  a 2x4 grid per global cell, showing the spatial trajectory and activity rate map
  or indicating if the cell was not tracked in that session.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from core.data_loader import list_sessions, load_session, mapping_rows_for_columns
from core.place_fields import gaussian_rate_maps_matlab


DEFAULT_SESSION_FILE = "2026_06_18_merged.mat"
DEFAULT_OUTPUT_DIR = Path("results/HabL_cell_evolution")
BIN_SIZE_CM = 1.0
SIGMA_CM = 4.0
ARENA_RADIUS_CM = 25.0


def precalculate_subsessions(data):
    """Precompute trajectories and continuous rate maps for all cells in all 4 subsessions."""
    mapping = np.asarray(data["mapping"], dtype=float)
    subsession_cache = []

    for s_idx, sub in enumerate(data["subsessions"]):
        x_raw = np.asarray(sub["x"], dtype=float)
        y_raw = np.asarray(sub["y"], dtype=float)
        t_raw = np.asarray(sub["t"], dtype=float)
        raw_S = np.asarray(sub["S"], dtype=float)

        valid = (
            np.isfinite(x_raw)
            & np.isfinite(y_raw)
            & np.isfinite(t_raw)
            & np.isfinite(raw_S).all(axis=1)
        )
        x = x_raw[valid]
        y = y_raw[valid]
        t = t_raw[valid]
        S_positive = np.maximum(raw_S[valid], 0.0).astype(np.float32)

        rates, edges, meta = gaussian_rate_maps_matlab(
            x,
            y,
            t,
            S_positive,
            bin_size=BIN_SIZE_CM,
            sigma=SIGMA_CM,
            arena_radius=ARENA_RADIUS_CM,
        )

        # Build map from global row -> local column index
        mapping_rows = mapping_rows_for_columns(mapping, s_idx, S_positive.shape[1])
        row_to_col = {int(row): col for col, row in enumerate(mapping_rows)}

        subsession_cache.append({
            "x": x,
            "y": y,
            "t": t,
            "S": S_positive,
            "rates": rates,
            "edges": edges,
            "sample_rate_hz": meta["sample_rate_hz"],
            "row_to_col": row_to_col,
        })

    return subsession_cache


def style_axis(ax):
    ax.set_aspect("equal")
    ax.set_xlim(-ARENA_RADIUS_CM, ARENA_RADIUS_CM)
    ax.set_ylim(-ARENA_RADIUS_CM, ARENA_RADIUS_CM)
    ax.set_xticks((-25, 0, 25))
    ax.set_yticks((-25, 0, 25))
    ax.tick_params(labelsize=6)


def draw_cell_evolution_figure(animal, global_cell, mapping_row, sub_cache):
    """Draw a 2x4 figure showing Trajectory (Row 1) and Rate Map (Row 2) across OF1..OF4."""
    sessions_present = [s_idx + 1 for s_idx in range(4) if np.isfinite(mapping_row[s_idx])]
    n_sessions = len(sessions_present)

    fig, axes = plt.subplots(2, 4, figsize=(12.8, 6.2), squeeze=False)

    cell_stats = {
        "animal": animal,
        "global_cell": int(global_cell),
        "n_sessions": n_sessions,
        "sessions_present": ",".join(f"OF{s}" for s in sessions_present),
    }

    for s_idx in range(4):
        ax_traj = axes[0, s_idx]
        ax_map = axes[1, s_idx]
        style_axis(ax_traj)
        style_axis(ax_map)

        local_id_val = mapping_row[s_idx]
        if np.isnan(local_id_val):
            # Not tracked in this subsession
            circle1 = plt.Circle((0, 0), ARENA_RADIUS_CM, color="0.85", fill=False, linestyle="--", linewidth=0.8)
            ax_traj.add_artist(circle1)
            ax_traj.text(0, 0, f"Not tracked in OF{s_idx + 1}", ha="center", va="center", color="0.55", fontsize=8)
            ax_traj.set_title(f"OF{s_idx + 1}: Not tracked", fontsize=8, color="0.45")

            circle2 = plt.Circle((0, 0), ARENA_RADIUS_CM, color="0.85", fill=False, linestyle="--", linewidth=0.8)
            ax_map.add_artist(circle2)
            ax_map.text(0, 0, f"Not tracked in OF{s_idx + 1}", ha="center", va="center", color="0.55", fontsize=8)
            ax_map.set_title(f"OF{s_idx + 1}: Not tracked", fontsize=8, color="0.45")

            cell_stats[f"OF{s_idx + 1}_local_id"] = np.nan
            cell_stats[f"OF{s_idx + 1}_S_gt_0_samples"] = 0
            cell_stats[f"OF{s_idx + 1}_peak_rate"] = np.nan
            continue

        local_id = int(local_id_val)
        cache = sub_cache[s_idx]
        local_col = cache["row_to_col"][int(global_cell)]

        x = cache["x"]
        y = cache["y"]
        sig = cache["S"][:, local_col]
        rate_map = cache["rates"][local_col]

        positive_mask = sig > 0
        n_positive = int(positive_mask.sum())
        peak_rate = float(np.nanmax(rate_map)) if not np.isnan(rate_map).all() else 0.0

        cell_stats[f"OF{s_idx + 1}_local_id"] = local_id
        cell_stats[f"OF{s_idx + 1}_S_gt_0_samples"] = n_positive
        cell_stats[f"OF{s_idx + 1}_peak_rate"] = peak_rate

        # 1. Trajectory with S > 0 scatter
        ax_traj.plot(x, y, color="0.82", linewidth=0.3, rasterized=True, zorder=1)
        if n_positive > 0:
            sc = ax_traj.scatter(
                x[positive_mask],
                y[positive_mask],
                c=sig[positive_mask],
                cmap="turbo",
                vmin=0,
                vmax=float(np.max(sig[positive_mask])),
                s=4.5,
                alpha=0.85,
                linewidths=0,
                rasterized=True,
                zorder=2,
            )
            cb_traj = plt.colorbar(sc, ax=ax_traj, fraction=0.046, pad=0.03)
            cb_traj.ax.tick_params(labelsize=6)
            cb_traj.set_label("S", fontsize=6)
        ax_traj.set_title(
            f"OF{s_idx + 1} (local {local_id})\n{n_positive} S > 0 samples",
            fontsize=8,
        )

        # 2. Activity rate map ("All recordings")
        im = ax_map.imshow(
            rate_map.T,
            origin="lower",
            extent=[-ARENA_RADIUS_CM, ARENA_RADIUS_CM, -ARENA_RADIUS_CM, ARENA_RADIUS_CM],
            cmap="turbo",
            vmin=0,
            vmax=peak_rate if peak_rate > 0 else 1.0,
            interpolation="nearest",
            rasterized=True,
        )
        cb_map = plt.colorbar(im, ax=ax_map, fraction=0.046, pad=0.03)
        cb_map.ax.tick_params(labelsize=6)
        cb_map.set_label("act/s", fontsize=6)
        ax_map.set_title(
            f"All recordings\nmax: {peak_rate:.2f} act/s",
            fontsize=8,
        )

    fig.suptitle(
        f"{animal} | Global cell {global_cell:04d} | Tracked in {n_sessions}/4 HabL sessions ({cell_stats['sessions_present']})",
        fontsize=11,
        fontweight="semibold",
        y=0.985,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94), w_pad=0.8, h_pad=0.8)
    return fig, cell_stats


_worker_sub_cache = None


def _init_worker(sub_cache):
    global _worker_sub_cache
    _worker_sub_cache = sub_cache


def _render_single_cell(args):
    animal, global_cell, mapping_row, animal_dir, overwrite = args
    sessions_present = [s_idx + 1 for s_idx in range(4) if np.isfinite(mapping_row[s_idx])]
    filename = f"{animal}_global_{global_cell:04d}_{len(sessions_present)}sessions.png"
    fig_path = animal_dir / filename

    if fig_path.exists() and not overwrite:
        stats = {
            "animal": animal,
            "global_cell": int(global_cell),
            "n_sessions": len(sessions_present),
            "sessions_present": ",".join(f"OF{s}" for s in sessions_present),
            "png": str(fig_path),
        }
        for s_idx in range(4):
            local_id_val = mapping_row[s_idx]
            if np.isnan(local_id_val):
                stats[f"OF{s_idx + 1}_local_id"] = np.nan
                stats[f"OF{s_idx + 1}_S_gt_0_samples"] = 0
                stats[f"OF{s_idx + 1}_peak_rate"] = np.nan
            else:
                cache = _worker_sub_cache[s_idx]
                local_col = cache["row_to_col"][int(global_cell)]
                sig = cache["S"][:, local_col]
                rate_map = cache["rates"][local_col]
                stats[f"OF{s_idx + 1}_local_id"] = int(local_id_val)
                stats[f"OF{s_idx + 1}_S_gt_0_samples"] = int((sig > 0).sum())
                stats[f"OF{s_idx + 1}_peak_rate"] = (
                    float(np.nanmax(rate_map)) if not np.isnan(rate_map).all() else 0.0
                )
        return stats

    fig, stats = draw_cell_evolution_figure(animal, global_cell, mapping_row, _worker_sub_cache)
    fig.savefig(fig_path, dpi=130, facecolor="white")
    plt.close(fig)
    stats["png"] = str(fig_path)
    return stats


def process_animal(animal, session_file, output_base, min_sessions=2, limit=None, max_workers=6, overwrite=False):
    """Process all eligible cells for a given animal with parallel rendering."""
    import concurrent.futures

    animal_dir = output_base / animal
    animal_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {animal} {session_file}...", flush=True)
    data = load_session(animal, session_file)
    mapping = np.asarray(data["mapping"], dtype=float)

    # Precalculate subsessions
    t0 = time.time()
    sub_cache = precalculate_subsessions(data)
    print(f"  Precomputed all 4 subsession rate maps in {time.time() - t0:.2f}s", flush=True)

    # Identify eligible global cells
    sessions_count = np.sum(np.isfinite(mapping), axis=1)
    eligible_indices = np.flatnonzero(sessions_count >= min_sessions)
    if limit is not None:
        eligible_indices = eligible_indices[:limit]
    total_cells = len(eligible_indices)
    print(f"  Found {total_cells} global cells present in >={min_sessions} sessions (limit={limit})", flush=True)

    tasks = [
        (animal, int(gc), mapping[gc], animal_dir, overwrite)
        for gc in eligible_indices
    ]

    rows = []
    t_start = time.time()
    if max_workers and max_workers > 1 and total_cells > 10:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=max_workers,
            initializer=_init_worker,
            initargs=(sub_cache,),
        ) as executor:
            for count, stats in enumerate(executor.map(_render_single_cell, tasks, chunksize=8), start=1):
                rows.append(stats)
                if count % 50 == 0 or count == total_cells:
                    elapsed = time.time() - t_start
                    speed = count / elapsed if elapsed > 0 else 0
                    print(f"  [{count}/{total_cells}] processed in {elapsed:.1f}s ({speed:.1f} cells/s)", flush=True)
    else:
        _init_worker(sub_cache)
        for count, task in enumerate(tasks, start=1):
            stats = _render_single_cell(task)
            rows.append(stats)
            if count % 50 == 0 or count == total_cells:
                elapsed = time.time() - t_start
                speed = count / elapsed if elapsed > 0 else 0
                print(f"  [{count}/{total_cells}] processed in {elapsed:.1f}s ({speed:.1f} cells/s)", flush=True)

    return rows


def generate_contact_sheets(animal_rows, output_base, top_n=24, per_page=4):
    """Create summary contact sheets for top-tracked and active cells."""
    from PIL import Image

    contact_dir = output_base / "contact_sheets"
    contact_dir.mkdir(parents=True, exist_ok=True)

    # Sort by sessions tracked (descending) and max peak rate (descending)
    def cell_score(r):
        peaks = [
            r.get(f"OF{i}_peak_rate", 0) for i in range(1, 5)
            if np.isfinite(r.get(f"OF{i}_peak_rate", 0))
        ]
        max_p = max(peaks) if peaks else 0.0
        return (r["n_sessions"], max_p)

    sorted_rows = sorted(animal_rows, key=cell_score, reverse=True)
    selected = [r for r in sorted_rows if Path(r["png"]).exists()][:top_n]

    if not selected:
        return

    pages = int(np.ceil(len(selected) / per_page))
    for page_idx in range(pages):
        chunk = selected[page_idx * per_page : (page_idx + 1) * per_page]
        images = [Image.open(r["png"]) for r in chunk]

        # Stack vertically
        max_w = max(img.width for img in images)
        total_h = sum(img.height for img in images)
        contact_img = Image.new("RGB", (max_w, total_h), color=(255, 255, 255))

        y_offset = 0
        for img in images:
            contact_img.paste(img, (0, y_offset))
            y_offset += img.height

        contact_path = contact_dir / f"contact_sheet_{page_idx + 1:02d}.png"
        contact_img.save(contact_path, quality=90)
        for img in images:
            img.close()

    print(f"Generated {pages} contact sheets in {contact_dir}", flush=True)


def animals_with_session(session_file):
    return sorted({
        session["animal"]
        for session in list_sessions()
        if session["session_file"] == session_file
    })


def main():
    parser = argparse.ArgumentParser(description="Generate HabL cell evolution figures.")
    parser.add_argument(
        "--animals",
        nargs="+",
        default=None,
        help="Animal IDs. If omitted, use every animal containing --session-file.",
    )
    parser.add_argument("--session-file", default=DEFAULT_SESSION_FILE)
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--min-sessions", type=int, default=2)
    parser.add_argument("--limit", type=int, default=None, help="Limit cells per animal for testing")
    parser.add_argument("--max-workers", type=int, default=6)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    animals = args.animals or animals_with_session(args.session_file)
    if not animals:
        raise SystemExit(
            f"No se encontró {args.session_file} para ningún animal en la carpeta de datos."
        )

    output_base = Path(args.output_dir)
    output_base.mkdir(parents=True, exist_ok=True)

    all_rows = []
    for animal in animals:
        print(f"\n==================== Processing {animal} ====================", flush=True)
        animal_rows = process_animal(
            animal,
            args.session_file,
            output_base,
            min_sessions=args.min_sessions,
            limit=args.limit,
            max_workers=args.max_workers,
            overwrite=args.overwrite,
        )
        all_rows.extend(animal_rows)

    # Save summary CSV
    summary_csv = output_base / "evolution_summary.csv"
    if all_rows:
        fieldnames = list(all_rows[0].keys())
        with summary_csv.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_rows)
        print(f"\nSaved summary of {len(all_rows)} cells to {summary_csv}", flush=True)

        generate_contact_sheets(all_rows, output_base, top_n=24, per_page=4)

    # Print summary breakdown
    print("\nSummary of cells generated:")
    for animal in animals:
        a_rows = [r for r in all_rows if r["animal"] == animal]
        in_4 = sum(1 for r in a_rows if r["n_sessions"] == 4)
        in_3 = sum(1 for r in a_rows if r["n_sessions"] == 3)
        in_2 = sum(1 for r in a_rows if r["n_sessions"] == 2)
        print(f"  {animal}: {len(a_rows)} cells total (in 4 sessions: {in_4}, in 3: {in_3}, in 2: {in_2})")


if __name__ == "__main__":
    main()
