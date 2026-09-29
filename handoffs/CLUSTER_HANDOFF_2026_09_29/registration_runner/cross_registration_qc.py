"""Centroid and field-of-view QC for a selected MiniAn registration run."""

import itertools
import json
import re
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd


def _session_columns(mappings):
    columns = [
        column for column in mappings.columns
        if isinstance(column, tuple) and len(column) > 1 and column[0] == "session"
    ]
    if not columns:
        raise ValueError("Mapping has no session columns")
    return columns


def _normalise_image(image):
    values = np.asarray(image, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return np.zeros(values.shape, dtype=float)
    lo, hi = np.percentile(finite, [2, 99])
    if not np.isfinite(hi - lo) or hi <= lo:
        return np.zeros(values.shape, dtype=float)
    scaled = np.clip((values - lo) / (hi - lo), 0, 1)
    scaled[~np.isfinite(scaled)] = 0
    return scaled


def _image_extent(widths, heights):
    return (
        float(widths.min()) - 0.5,
        float(widths.max()) + 0.5,
        float(heights.max()) + 0.5,
        float(heights.min()) - 0.5,
    )


def _draw_mask(ax, widths, heights, mask, color, linestyle, linewidth):
    mask = np.asarray(mask, dtype=bool)
    if mask.any() and not mask.all():
        ax.contour(
            widths, heights, mask.astype(float), levels=[0.5],
            colors=[color], linestyles=[linestyle], linewidths=linewidth,
        )


def _short_day(label):
    day = label.split("/")[0]
    day = re.split(r"_?20\d{2}_\d{2}_\d{2}", day)[0]
    return day.strip(" _") or label.split("/")[0]


def _slug(label):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", label).strip("_")


def _centroid_lookup(cents):
    lookup = {}
    for row in cents[["session", "unit_id", "width", "height"]].itertuples(index=False):
        key = (str(row.session), int(row.unit_id))
        if key in lookup:
            raise ValueError("Duplicate centroid for {} unit_id {}".format(*key))
        lookup[key] = (float(row.width), float(row.height))
    return lookup


def _pair_matches(mappings, session_a, session_b, centroid_lookup):
    values_a = mappings[("session", session_a)]
    values_b = mappings[("session", session_b)]
    matched = values_a.notna() & values_b.notna()
    rows = []
    for global_id in mappings.index[matched]:
        unit_a = int(values_a.loc[global_id])
        unit_b = int(values_b.loc[global_id])
        x_a, y_a = centroid_lookup[(session_a, unit_a)]
        x_b, y_b = centroid_lookup[(session_b, unit_b)]
        rows.append({
            "global_cell_id": int(global_id),
            "unit_id_a": unit_a,
            "unit_id_b": unit_b,
            "x_a": x_a,
            "y_a": y_a,
            "x_b": x_b,
            "y_b": y_b,
            "centroid_distance_pixels": float(np.hypot(x_a - x_b, y_a - y_b)),
        })
    return rows


def _save_overview(
    folder, groups, sessions, display, cents, mappings, shifted,
    common_mask, widths, heights, output_dir,
):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    n_rows = len(groups)
    n_cols = max(len(items) for items in groups.values())
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(4.1 * n_cols, 4.1 * n_rows), squeeze=False
    )
    multi_day = np.zeros(len(mappings), dtype=int)
    for group_sessions in groups.values():
        group_columns = [("session", session) for session in group_sessions]
        multi_day += mappings[group_columns].notna().any(axis=1).to_numpy(dtype=int)
    shared_ids = set(mappings.index[multi_day >= 2])
    extent = _image_extent(widths, heights)

    for row_number, (day, group_sessions) in enumerate(groups.items()):
        for column_number, session in enumerate(group_sessions):
            ax = axes[row_number, column_number]
            image = shifted[session]
            ax.imshow(
                np.ma.masked_invalid(image), cmap="gray", extent=extent,
                origin="upper", interpolation="nearest",
            )
            _draw_mask(ax, widths, heights, common_mask, "#ffcf40", "-", 1.1)
            local = cents.loc[cents["session"] == session]
            ax.scatter(
                local["width"], local["height"], s=2, c="#55d6ef",
                alpha=0.45, linewidths=0,
            )
            shared_units = [
                int(value) for global_id, value in mappings[("session", session)].items()
                if global_id in shared_ids and pd.notna(value)
            ]
            marked = local.loc[local["unit_id"].isin(shared_units)]
            ax.scatter(
                marked["width"], marked["height"], s=10,
                facecolors="none", edgecolors="#ff8a4a", linewidths=0.6,
            )
            ax.set_title(
                "{} / {}\n{} centroids; {} linked across days".format(
                    _short_day(display[session]), display[session].split("/")[-1],
                    len(local), len(marked),
                ),
                fontsize=8,
            )
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
            ax.set_aspect("equal", adjustable="box")
            ax.set_xticks([])
            ax.set_yticks([])
        for column_number in range(len(group_sessions), n_cols):
            axes[row_number, column_number].axis("off")

    handles = [
        Line2D([0], [0], color="#ffcf40", lw=1.2, label="common FOV (all selected sessions)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#55d6ef",
               markersize=5, label="usable centroid"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="none",
               markeredgecolor="#ff8a4a", markersize=6, label="cell linked across days"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=9)
    fig.suptitle("Aligned fields of view and centroids: {}".format(folder), fontsize=14)
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    path = output_dir / "centroid_fov_overview.png"
    fig.savefig(str(path), dpi=150)
    plt.close(fig)
    return path


def _save_pair_panels(
    day_a, day_b, sessions_a, sessions_b, display, cents, pair_details,
    pair_stats, shifted, common_mask, widths, heights, output_dir, distance_limit,
):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    same_day = day_a == day_b
    n_rows, n_cols = len(sessions_a), len(sessions_b)
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(4.1 * n_cols, 4.1 * n_rows), squeeze=False
    )
    extent = _image_extent(widths, heights)
    for row_number, session_a in enumerate(sessions_a):
        for column_number, session_b in enumerate(sessions_b):
            ax = axes[row_number, column_number]
            if same_day and column_number <= row_number:
                ax.axis("off")
                continue
            image_a = shifted[session_a]
            image_b = shifted[session_b]
            pair_mask = np.isfinite(image_a) & np.isfinite(image_b)
            rgb = np.zeros(image_a.shape + (3,), dtype=float)
            rgb[:, :, 0] = _normalise_image(image_a)
            rgb[:, :, 1] = _normalise_image(image_b)
            rgb[~pair_mask] *= 0.2
            ax.imshow(rgb, extent=extent, origin="upper", interpolation="nearest")
            _draw_mask(ax, widths, heights, pair_mask, "#a7f08d", "-", 0.8)
            _draw_mask(ax, widths, heights, common_mask, "#ffffff", "--", 1.0)

            cells_a = cents.loc[cents["session"] == session_a]
            cells_b = cents.loc[cents["session"] == session_b]
            ax.scatter(
                cells_a["width"], cells_a["height"], s=2, c="#ff66aa",
                alpha=0.25, linewidths=0,
            )
            ax.scatter(
                cells_b["width"], cells_b["height"], s=2, c="#56d9ff",
                alpha=0.25, linewidths=0,
            )
            details = pair_details[(session_a, session_b)]
            for item in details:
                color = "#ff9b54" if item["centroid_distance_pixels"] > distance_limit else "#ffffff"
                ax.plot(
                    [item["x_a"], item["x_b"]], [item["y_a"], item["y_b"]],
                    color=color, alpha=0.52, linewidth=0.55,
                )
            if details:
                ax.scatter(
                    [item["x_a"] for item in details],
                    [item["y_a"] for item in details],
                    s=10, c="#ff66aa", edgecolors="black", linewidths=0.2,
                )
                ax.scatter(
                    [item["x_b"] for item in details],
                    [item["y_b"] for item in details],
                    s=12, c="#56d9ff", marker="^",
                    edgecolors="black", linewidths=0.2,
                )
            stats = pair_stats[(session_a, session_b)]
            ax.set_title(
                "{} vs {}\n{} linked | {:.1f}% / {:.1f}% | {} > {}px".format(
                    display[session_a].split("/")[-1],
                    display[session_b].split("/")[-1],
                    stats["n_matched"], stats["pct_a"], stats["pct_b"],
                    stats["n_over_distance_limit"], distance_limit,
                ),
                fontsize=7,
            )
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
            ax.set_aspect("equal", adjustable="box")
            ax.set_xticks([])
            ax.set_yticks([])

    handles = [
        Line2D([0], [0], color="#a7f08d", lw=1, label="pair overlap area"),
        Line2D([0], [0], color="black", ls="--", lw=1, label="common area across all sessions"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#ff66aa",
               markersize=5, label="session A centroid"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor="#56d9ff",
               markersize=5, label="session B centroid"),
        Line2D([0], [0], color="#ff9b54", lw=1, label="matched pair over distance limit"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8)
    fig.suptitle(
        "Aligned areas and matched centroids: {} vs {}".format(day_a, day_b),
        fontsize=13,
    )
    fig.tight_layout(rect=(0, 0.07, 1, 0.96))
    path = output_dir / (
        "centroid_matches_{}__{}.png".format(_slug(day_a), _slug(day_b))
    )
    fig.savefig(str(path), dpi=160)
    plt.close(fig)
    return path


def create_registration_qc(
    mappings, cents, shiftds, common_window_2d, session_name_map,
    output_dir, logger, distance_limit, make_figures,
):
    """Save matched-cell percentages, FOV areas, and centroid overlays."""
    output_dir = Path(output_dir)
    sessions = [column[1] for column in _session_columns(mappings)]
    display = {
        session: session_name_map.get(session, session) for session in sessions
    }
    groups = OrderedDict()
    for session in sessions:
        day = display[session].split("/")[0]
        groups.setdefault(day, []).append(session)
    centroid_lookup = _centroid_lookup(cents)
    n_centroids = cents.groupby("session")["unit_id"].nunique().to_dict()

    shifted_array = shiftds["temps_shifted"]
    shifted = {
        session: np.asarray(
            shifted_array.sel(session=session).transpose("height", "width").values,
            dtype=float,
        )
        for session in sessions
    }
    widths = np.asarray(shifted_array.coords["width"].values, dtype=float)
    heights = np.asarray(shifted_array.coords["height"].values, dtype=float)
    common_mask = np.asarray(
        common_window_2d.transpose("height", "width").compute().values, dtype=bool
    )
    n_total_pixels = int(common_mask.size)
    n_common_pixels = int(common_mask.sum())

    pair_details = {}
    pair_stats = {}
    summary_rows = []
    detail_rows = []
    for session_a, session_b in itertools.combinations(sessions, 2):
        details = _pair_matches(mappings, session_a, session_b, centroid_lookup)
        pair_details[(session_a, session_b)] = details
        distances = np.asarray(
            [item["centroid_distance_pixels"] for item in details], dtype=float
        )
        pair_fov_pixels = int(
            (np.isfinite(shifted[session_a]) & np.isfinite(shifted[session_b])).sum()
        )
        n_a = int(n_centroids.get(session_a, 0))
        n_b = int(n_centroids.get(session_b, 0))
        n_matched = len(details)
        stats = {
            "day_a": display[session_a].split("/")[0],
            "session_a": display[session_a],
            "day_b": display[session_b].split("/")[0],
            "session_b": display[session_b],
            "n_centroids_a": n_a,
            "n_centroids_b": n_b,
            "n_matched": n_matched,
            "pct_a": 100.0 * n_matched / n_a if n_a else float("nan"),
            "pct_b": 100.0 * n_matched / n_b if n_b else float("nan"),
            "median_distance_pixels": float(np.median(distances)) if distances.size else float("nan"),
            "p95_distance_pixels": float(np.percentile(distances, 95)) if distances.size else float("nan"),
            "n_over_distance_limit": int((distances > distance_limit).sum()),
            "pair_fov_pixels": pair_fov_pixels,
            "pair_fov_pct_full": 100.0 * pair_fov_pixels / n_total_pixels,
            "all_selected_fov_pixels": n_common_pixels,
            "all_selected_fov_pct_full": 100.0 * n_common_pixels / n_total_pixels,
        }
        pair_stats[(session_a, session_b)] = stats
        summary_rows.append(stats)
        for item in details:
            detail_rows.append(dict(stats, **item))
    pd.DataFrame(summary_rows).to_csv(
        output_dir / "pairwise_match_summary.csv", index=False
    )
    pd.DataFrame(detail_rows).to_csv(
        output_dir / "centroid_match_pairs.csv.gz",
        index=False, compression="gzip",
    )

    day_presence = {}
    for day, group_sessions in groups.items():
        group_columns = [("session", session) for session in group_sessions]
        day_presence[day] = mappings[group_columns].notna().any(axis=1)
    day_rows = []
    for day_a, day_b in itertools.combinations(groups, 2):
        present_a, present_b = day_presence[day_a], day_presence[day_b]
        n_a, n_b = int(present_a.sum()), int(present_b.sum())
        shared = int((present_a & present_b).sum())
        day_rows.append({
            "day_a": day_a, "day_b": day_b,
            "n_global_ids_a": n_a, "n_global_ids_b": n_b,
            "n_shared_global_ids": shared,
            "pct_a": 100.0 * shared / n_a if n_a else float("nan"),
            "pct_b": 100.0 * shared / n_b if n_b else float("nan"),
        })
    pd.DataFrame(day_rows).to_csv(output_dir / "day_match_summary.csv", index=False)
    all_days_shared = int(np.logical_and.reduce(
        [presence.to_numpy(dtype=bool) for presence in day_presence.values()]
    ).sum())
    scope = {
        "sessions": [display[session] for session in sessions],
        "n_sessions": len(sessions),
        "day_groups": list(groups),
        "n_all_selected_fov_pixels": n_common_pixels,
        "all_selected_fov_pct_full": 100.0 * n_common_pixels / n_total_pixels,
        "n_global_ids_present_in_every_day": all_days_shared,
        "n_global_ids_present_in_every_acquisition": int(
            mappings[[("session", session) for session in sessions]].notna().all(axis=1).sum()
        ),
    }
    (output_dir / "registration_scope.json").write_text(
        json.dumps(scope, indent=2), encoding="utf-8"
    )
    logger.info("Common FOV: %d / %d pixels (%.1f%%)", n_common_pixels,
                n_total_pixels, scope["all_selected_fov_pct_full"])
    logger.info("Global IDs seen in every selected day: %d", all_days_shared)
    logger.info("Saved pairwise and day-level match percentages to %s", output_dir)

    figure_paths = []
    if make_figures:
        overview_path = _save_overview(
            output_dir.name, groups, sessions, display, cents, mappings,
            shifted, common_mask, widths, heights, output_dir,
        )
        figure_paths.append(overview_path)
        for day_a, day_b in itertools.combinations_with_replacement(groups, 2):
            path = _save_pair_panels(
                day_a, day_b, groups[day_a], groups[day_b],
                display, cents, pair_details, pair_stats, shifted,
                common_mask, widths, heights, output_dir, distance_limit,
            )
            figure_paths.append(path)
        for path in figure_paths:
            logger.info("Saved centroid/FOV QC figure: %s", path)
    return scope, summary_rows, day_rows, figure_paths

