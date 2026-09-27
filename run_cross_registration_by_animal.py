#!/usr/bin/env python3
"""
Cross-session cell registration for MiniAn outputs from one animal.

Expected folder structure:

    ANIMAL_ROOT/
        day_or_protocol_1/acquisition_1/My_V4_Miniscope/minian/
        day_or_protocol_1/acquisition_2/My_V4_Miniscope/minian/
        day_or_protocol_2/acquisition_1/My_V4_Miniscope/minian/
        ...

Pass one animal folder with --animal-root. The script recursively finds each
My_V4_Miniscope folder under it and uses the relative day/acquisition path as
the unique session label. This prevents repeated names such as A_AV1 from
colliding across different days. Every discovered My_V4_Miniscope folder must
contain an existing minian folder. If any is missing, the script stops.

Outputs are saved in a new timestamped folder under ANIMAL_ROOT by default:

    ANIMAL_ROOT/cross_registration_YYYYMMDD_HHMMSS/...

The mapping index is a new, consecutive ``global_cell_id`` (1..N).  Values in
the session columns are the original MiniAn ``unit_id`` labels and are not
required to be consecutive.
"""

# =============================================================================
# USER PARAMETERS
# =============================================================================

from pathlib import Path

# Folder names expected inside each session.
MINISCOPE_FOLDER_NAME = "My_V4_Miniscope"
MINIAN_FOLDER_NAME = "minian"

# Main cross-registration parameter.
# Maximum centroid distance, in pixels, to consider cells from different
# sessions as candidates for the same cell.
PARAM_DIST_PIXELS = 5

# Output names saved in the new output directory.
LOG_FILE_NAME = "cross_registration_log.log"
MAPPINGS_FILE_NAME = "mappings.pkl"
MAPPINGS_CSV_FILE_NAME = "mappings.csv"
MATLAB_MAPPINGS_CSV_FILE_NAME = "mappings_matlab.csv"
UNIT_ID_LOOKUP_CSV_FILE_NAME = "unit_id_to_C_column.csv"
CENTROIDS_FILE_NAME = "cents.pkl"
SHIFTS_FILE_NAME = "shiftds.nc"
ALIGNMENT_FIGURE_NAME = "cross_registration_alignment.png"
CONTOURS_FILE_NAME = "cell_contours.csv.gz"
CONTOUR_FIGURE_NAME = "cross_registration_contours.pdf"

# Only remove known outputs inside the selected output directory.
OVERWRITE_OUTPUTS = True

# Static image options.
FIGURE_DPI = 150
FIGURE_ROW_HEIGHT = 3.0
FIGURE_COL_WIDTH = 4.0

# Compact matched-cell QC.  The CSV contains only vector contour coordinates
# (plus identifiers/centroids) and is gzip-compressed.  The PDF places many
# local contour overlays on each page instead of saving one full-frame image
# per cell.
SAVE_CONTOUR_QC = True
CONTOUR_LEVEL_FRACTION = 0.20
CONTOUR_KEEP_LARGEST_ONLY = True
CONTOUR_MIN_SESSIONS = 2
CONTOUR_PANELS_PER_ROW = 5
CONTOUR_ROWS_PER_PAGE = 6
CONTOUR_PADDING_PIXELS = 4

# =============================================================================
# IMPORTS
# =============================================================================

import argparse
import logging
import os
import sys
import time
import traceback
from collections import Counter
from contextlib import contextmanager

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.lines import Line2D
    HAS_MATPLOTLIB = True
except ImportError:
    plt = None
    PdfPages = None
    Line2D = None
    HAS_MATPLOTLIB = False
import numpy as np
import pandas as pd
import xarray as xr
from dask.diagnostics import ProgressBar
from skimage.measure import find_contours

from minian.cross_registration import (
    calculate_centroids,
    calculate_centroid_distance,
    calculate_mapping,
    fill_mapping,
    group_by_session,
    resolve_mapping,
)
from minian.motion_correction import apply_transform, estimate_motion
from minian.utilities import open_minian


# =============================================================================
# LOGGING AND TIMING HELPERS
# =============================================================================

def setup_logger(log_file):
    """Create one logger that writes to terminal and to the output folder."""
    log_file = Path(log_file)
    log_file.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger("cross_registration")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False

    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(fmt)
    stream_handler.setLevel(logging.INFO)
    logger.addHandler(stream_handler)

    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setFormatter(fmt)
    file_handler.setLevel(logging.INFO)
    logger.addHandler(file_handler)

    return logger


@contextmanager
def timed_step(logger, message):
    """Log wall time for a processing step."""
    logger.info(message)
    t0 = time.time()
    try:
        yield
    finally:
        logger.info("%s done in %.1f s", message, time.time() - t0)


# =============================================================================
# FOLDER DISCOVERY AND VALIDATION
# =============================================================================

def find_session_miniscope_folders(dpath):
    """Find My_V4_Miniscope folders, pruning each match from the walk."""
    dpath = Path(dpath)
    found = []
    for current, dirnames, _ in os.walk(dpath):
        current_path = Path(current)
        if current_path.name == MINISCOPE_FOLDER_NAME:
            found.append(current_path)
            # Do not descend into minian/Zarr contents.
            dirnames[:] = []
    return sorted(found)


def validate_minian_folders(session_folders):
    """Stop if any discovered session lacks a minian output folder."""
    missing = []
    for session_path in session_folders:
        minian_path = session_path / MINIAN_FOLDER_NAME
        if not minian_path.is_dir():
            missing.append(minian_path)

    if missing:
        msg = [
            "The following sessions have My_V4_Miniscope folders but no minian folder:",
            "",
        ]
        msg.extend([os.fspath(p) for p in missing])
        raise FileNotFoundError("\n".join(msg))


def remove_previous_outputs(dpath, logger):
    """Remove previous cross-registration outputs when overwriting is enabled."""
    if not OVERWRITE_OUTPUTS:
        return

    for name in [
        MAPPINGS_FILE_NAME,
        MAPPINGS_CSV_FILE_NAME,
        MATLAB_MAPPINGS_CSV_FILE_NAME,
        UNIT_ID_LOOKUP_CSV_FILE_NAME,
        CENTROIDS_FILE_NAME,
        SHIFTS_FILE_NAME,
        ALIGNMENT_FIGURE_NAME,
        CONTOURS_FILE_NAME,
        CONTOUR_FIGURE_NAME,
    ]:
        path = Path(dpath) / name
        if path.exists():
            path.unlink()
            logger.info("Removed previous output: %s", path)


# =============================================================================
# DATA LOADING
# =============================================================================

def open_all_sessions(session_folders, logger, animal_root):
    """Open all MiniAn datasets and concatenate them along session."""
    datasets = []
    original_cells = {}
    original_cell_ids = {}
    c_column_lookups = {}
    session_name_map = {}
    reference_spatial_coords = None

    for session_number, session_path in enumerate(session_folders):
        session_path = Path(session_path)
        # Include both the day/protocol folder and acquisition folder. Leaf
        # names such as A_AV1 repeat across days, so using only the leaf name
        # would create duplicate mapping columns.
        display_name = session_path.parent.relative_to(animal_root).as_posix()
        # MiniAn's resolve_mapping currently encodes graph nodes as
        # ``session-unit_id`` and later splits on "-".  Safe internal labels
        # prevent folder names containing hyphens from corrupting that step.
        session_name = f"s{session_number:04d}"
        session_name_map[session_name] = display_name
        minian_path = session_path / MINIAN_FOLDER_NAME

        logger.info("Opening session: %s (internal label %s)", display_name, session_name)
        logger.info("MiniAn folder: %s", minian_path)

        # Loading the arrays separately avoids xr.merge aligning A and C onto
        # the union of their unit_id coordinates and padding missing cells with
        # NaN.  C is the authoritative filtered-cell set for downstream use.
        minian_vars = open_minian(os.fspath(minian_path), return_dict=True)

        for required_var in ["A", "C", "max_proj"]:
            if required_var not in minian_vars:
                raise KeyError(
                    f"Required variable '{required_var}' not found in {minian_path}"
                )
        A = minian_vars["A"]
        C = minian_vars["C"]
        max_proj = minian_vars["max_proj"]

        for dim in ["unit_id", "height", "width"]:
            if dim not in A.dims:
                raise ValueError(f"A in {minian_path} has no '{dim}' dimension")

        a_unit_ids = np.asarray(A.coords["unit_id"].values)
        if pd.isna(a_unit_ids).any():
            raise ValueError(f"NaN unit_id found in {minian_path}")
        if pd.Index(a_unit_ids).has_duplicates:
            duplicated = pd.Index(a_unit_ids)[
                pd.Index(a_unit_ids).duplicated()
            ].unique()
            raise ValueError(
                f"Duplicate unit_id values in {minian_path}: {duplicated.tolist()}"
            )
        try:
            a_unit_ids_float = a_unit_ids.astype(float)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"MiniAn cross-registration requires numeric unit_id values in {minian_path}"
            ) from exc
        if not np.all(np.isfinite(a_unit_ids_float)) or not np.allclose(
            a_unit_ids_float, np.round(a_unit_ids_float)
        ):
            raise ValueError(
                f"MiniAn cross-registration requires finite integer-like unit_id "
                f"values in {minian_path}"
            )

        if "unit_id" not in C.dims:
            raise ValueError(f"C in {minian_path} has no 'unit_id' dimension")
        c_unit_ids = np.asarray(C.coords["unit_id"].values)
        if pd.isna(c_unit_ids).any() or pd.Index(c_unit_ids).has_duplicates:
            raise ValueError(f"C has missing or duplicate unit_id values in {minian_path}")
        try:
            c_unit_ids_float = c_unit_ids.astype(float)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"C has non-numeric unit_id values in {minian_path}") from exc
        if not np.all(np.isfinite(c_unit_ids_float)) or not np.allclose(
            c_unit_ids_float, np.round(c_unit_ids_float)
        ):
            raise ValueError(
                f"C has non-integer-like unit_id values in {minian_path}"
            )
        a_unit_id_set = {int(value) for value in a_unit_ids_float}
        c_unit_id_set = {int(value) for value in c_unit_ids_float}
        only_c = sorted(c_unit_id_set - a_unit_id_set)
        if only_c:
            raise ValueError(
                f"C contains unit_id values without spatial footprints in A "
                f"in {minian_path}: {only_c[:10]}"
            )
        only_a = sorted(a_unit_id_set - c_unit_id_set)
        if only_a:
            logger.warning(
                "%s contains %d cells in A that are absent from filtered C. "
                "They will not be included in cross-registration. First IDs: %s",
                display_name,
                len(only_a),
                only_a[:10],
            )
        # MATLAB indices are 1-based.  This lookup is valid when act.C is made
        # from C.values without reordering its unit_id dimension.
        c_column_lookups[session_name] = {
            int(unit_id): column_number
            for column_number, unit_id in enumerate(c_unit_ids_float, start=1)
        }

        spatial_coords = {
            dim: np.asarray(A.coords[dim].values) for dim in ["height", "width"]
        }
        if reference_spatial_coords is None:
            reference_spatial_coords = spatial_coords
        else:
            for dim in ["height", "width"]:
                if not np.array_equal(spatial_coords[dim], reference_spatial_coords[dim]):
                    raise ValueError(
                        "All sessions must use the same spatial pixel grid. "
                        f"Coordinate '{dim}' differs in {minian_path}."
                    )

        n_cells = len(c_unit_ids)
        original_cells[display_name] = n_cells
        original_cell_ids[display_name] = c_unit_ids.copy()
        uid_min = int(np.min(c_unit_ids_float)) if n_cells else None
        uid_max = int(np.max(c_unit_ids_float)) if n_cells else None
        if n_cells:
            n_gaps = int(uid_max - uid_min + 1 - n_cells)
            logger.info(
                "Filtered C cells used in %s: %d; unit_id range %d..%d (%d gaps)",
                display_name,
                n_cells,
                uid_min,
                uid_max,
                n_gaps,
            )
            logger.info(
                "C columns in %s: %d (MATLAB lookup saved as 1-based indices)",
                display_name,
                len(c_unit_ids),
            )
        else:
            logger.info("Filtered C cells used in %s: 0", display_name)

        # Select A in exactly the C unit order.  This both excludes any
        # unfiltered A-only cells and makes the 1-based MATLAB lookup explicit.
        ds = xr.Dataset(
            {
                "A": A.sel(unit_id=c_unit_ids.tolist()),
                "max_proj": max_proj,
            }
        )
        ds = ds.expand_dims(session=[session_name])
        datasets.append(ds)

    minian_ds = xr.concat(datasets, dim="session", join="outer")

    logger.info("Sessions loaded: %s", list(minian_ds.coords["session"].values))
    logger.info("A sizes after concat: %s", dict(minian_ds["A"].sizes))

    return (
        minian_ds,
        original_cells,
        original_cell_ids,
        session_name_map,
        c_column_lookups,
    )


# =============================================================================
# ALIGNMENT AND STATIC QC FIGURE
# =============================================================================

def estimate_session_shifts(minian_ds):
    """Estimate translational shifts using max projections."""
    temps = minian_ds["max_proj"].rename("temps")
    shifts = estimate_motion(temps, dim="session").compute().rename("shifts")
    temps_shifted = apply_transform(temps, shifts).compute().rename("temps_shifted")
    shiftds = xr.merge([temps, shifts, temps_shifted])
    return shiftds


def make_common_window(shiftds):
    """Create the common field-of-view window used for centroid calculation."""
    # A true common FOV must be valid in every session.  Comparing with the
    # minimum NaN count is unsafe: if no pixel is shared by all sessions, that
    # expression silently selects a merely "best overlap" region.
    common_window_2d = shiftds["temps_shifted"].notnull().all("session")
    n_common = int(common_window_2d.sum().compute().item())
    if n_common == 0:
        raise ValueError(
            "The shifted sessions have no common field of view. "
            "Inspect the alignment figure and session shifts."
        )

    window, _ = xr.broadcast(common_window_2d, shiftds["temps_shifted"])
    return window, common_window_2d


def save_alignment_figure(
    shiftds,
    common_window_2d,
    figure_path,
    logger,
    session_name_map=None,
):
    """Save a static before/after alignment figure for sanity checking."""
    sessions = list(shiftds.coords["session"].values)
    n_sessions = len(sessions)

    fig_width = FIGURE_COL_WIDTH * 3
    fig_height = max(FIGURE_ROW_HEIGHT * n_sessions, 3.0)
    fig, axes = plt.subplots(
        n_sessions,
        3,
        figsize=(fig_width, fig_height),
        squeeze=False,
    )

    window_mask = common_window_2d.transpose("height", "width").values.astype(bool)

    for row, session in enumerate(sessions):
        display_session = (
            session_name_map.get(session, session) if session_name_map else session
        )
        original = (
            shiftds["temps"]
            .sel(session=session)
            .transpose("height", "width")
            .values
        )
        shifted = (
            shiftds["temps_shifted"]
            .sel(session=session)
            .transpose("height", "width")
            .values
        )
        shifted_common = np.array(shifted, copy=True)
        shifted_common[~window_mask] = np.nan

        images = [original, shifted, shifted_common]
        titles = [
            f"{display_session}\noriginal",
            f"{display_session}\nshifted",
            f"{display_session}\ncommon FOV",
        ]

        for col, (img, title) in enumerate(zip(images, titles)):
            ax = axes[row, col]
            ax.imshow(img, cmap="viridis")
            ax.set_title(title, fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])

    fig.tight_layout()
    figure_path = Path(figure_path)
    fig.savefig(figure_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved alignment QC figure: %s", figure_path)


# =============================================================================
# CROSS-SESSION REGISTRATION
# =============================================================================

def calculate_cross_session_mapping(minian_ds, shiftds, window, logger):
    """Calculate centroids, pairwise distances, and cross-session mappings."""
    with timed_step(logger, "Applying shifts to spatial footprints"):
        A_shifted = apply_transform(
            minian_ds["A"].chunk(dict(height=-1, width=-1)),
            shiftds["shifts"],
        )

    with timed_step(logger, "Calculating centroids"):
        cents = calculate_centroids(A_shifted, window)
        if cents.empty:
            raise ValueError("No usable cell centroids remain inside the common FOV")

    with timed_step(logger, "Calculating centroid distances"):
        # There are no additional identifying dimensions besides session in this
        # folder structure, so index_dim is empty.
        dist = calculate_centroid_distance(cents, index_dim=[])

    with timed_step(logger, f"Thresholding distances at {PARAM_DIST_PIXELS} pixels"):
        dist_filtered = dist[dist["variable", "distance"] < PARAM_DIST_PIXELS].copy()
        dist_filtered = group_by_session(dist_filtered)

    with timed_step(logger, "Generating mappings"):
        if dist_filtered.empty:
            logger.warning(
                "No cross-session centroid pairs were closer than %.2f pixels; "
                "all usable cells will be recorded as unmatched",
                PARAM_DIST_PIXELS,
            )
            sessions = list(minian_ds.coords["session"].values)
            rows = []
            for cent in cents[["session", "unit_id"]].itertuples(index=False):
                row = {("session", session): np.nan for session in sessions}
                row[("session", cent.session)] = cent.unit_id
                rows.append(row)
            mappings_meta_fill = pd.DataFrame.from_records(rows)
            mappings_meta_fill.columns = pd.MultiIndex.from_tuples(
                mappings_meta_fill.columns
            )
            mappings_meta_fill = group_by_session(mappings_meta_fill)
        else:
            mappings = calculate_mapping(dist_filtered)
            mappings_meta = resolve_mapping(mappings)
            mappings_meta_fill = fill_mapping(mappings_meta, cents)

            # A session with no usable centroid can otherwise disappear from
            # the mapping columns, making an "all sessions" count misleading.
            sessions = list(minian_ds.coords["session"].values)
            for session in sessions:
                col = ("session", session)
                if col not in mappings_meta_fill.columns:
                    mappings_meta_fill[col] = np.nan
            session_columns = [("session", session) for session in sessions]
            other_columns = [
                col for col in mappings_meta_fill.columns if col[0] != "session"
            ]
            mappings_meta_fill = mappings_meta_fill.reindex(
                columns=pd.MultiIndex.from_tuples(session_columns + other_columns)
            )
            mappings_meta_fill = group_by_session(mappings_meta_fill)

    return A_shifted, cents, mappings_meta_fill


# =============================================================================
# REPORTING AND SAVING
# =============================================================================

def get_session_columns(mappings):
    """Return mapping columns that contain original per-session unit IDs."""
    session_cols = [
        col
        for col in mappings.columns
        if isinstance(col, tuple) and len(col) > 0 and col[0] == "session"
    ]
    if not session_cols:
        raise ValueError(
            "Could not find session columns in mappings output. "
            f"Columns were: {list(mappings.columns)}"
        )
    return session_cols


def normalise_unit_id(value):
    """Convert MiniAn's float representation back to an integer label."""
    value_float = float(value)
    if not np.isfinite(value_float) or not value_float.is_integer():
        raise ValueError(f"Invalid mapped unit_id: {value!r}")
    return int(value_float)


def validate_mapping_integrity(mappings, cents):
    """Verify that every usable (session, unit_id) occurs exactly once."""
    session_cols = get_session_columns(mappings)
    expected = {
        (row.session, normalise_unit_id(row.unit_id))
        for row in cents[["session", "unit_id"]].itertuples(index=False)
    }
    observed_list = []
    for col in session_cols:
        session = col[1]
        observed_list.extend(
            (session, normalise_unit_id(uid))
            for uid in mappings[col].dropna().tolist()
        )

    observed_counts = Counter(observed_list)
    observed = set(observed_counts)
    duplicates = [cell for cell, count in observed_counts.items() if count != 1]
    missing = sorted(expected - observed, key=str)
    extra = sorted(observed - expected, key=str)
    if duplicates or missing or extra:
        raise ValueError(
            "Mapping integrity check failed. "
            f"duplicates={duplicates[:10]}, missing={missing[:10]}, extra={extra[:10]}"
        )


def add_mapping_qc_and_global_ids(mappings, cents, logger):
    """Add distance diagnostics and a deterministic consecutive global ID."""
    mappings = mappings.copy()
    session_cols = get_session_columns(mappings)
    cent_lookup = cents.set_index(["session", "unit_id"])[["height", "width"]]
    if cent_lookup.index.has_duplicates:
        raise ValueError("Centroid table has duplicate (session, unit_id) entries")

    n_sessions_list = []
    max_distance_list = []
    mean_distance_list = []
    mean_height_list = []
    mean_width_list = []

    for _, row in mappings.iterrows():
        coords = []
        for col in session_cols:
            if pd.isna(row[col]):
                continue
            session = col[1]
            unit_id = normalise_unit_id(row[col])
            try:
                cent = cent_lookup.loc[(session, unit_id)]
            except KeyError as exc:
                raise ValueError(
                    f"Mapping contains ({session}, {unit_id}) but no centroid exists"
                ) from exc
            coords.append((float(cent["height"]), float(cent["width"])))

        if not coords:
            raise ValueError("A mapping row contains no cells")
        coord_array = np.asarray(coords, dtype=float)
        pairwise = []
        for idx_a in range(len(coord_array)):
            for idx_b in range(idx_a + 1, len(coord_array)):
                pairwise.append(float(np.linalg.norm(coord_array[idx_a] - coord_array[idx_b])))

        n_sessions_list.append(len(coords))
        max_distance_list.append(max(pairwise) if pairwise else 0.0)
        mean_distance_list.append(float(np.mean(pairwise)) if pairwise else 0.0)
        mean_height_list.append(float(coord_array[:, 0].mean()))
        mean_width_list.append(float(coord_array[:, 1].mean()))

    mappings[("qc", "n_sessions")] = n_sessions_list
    mappings[("qc", "max_centroid_distance_pixels")] = max_distance_list
    mappings[("qc", "mean_centroid_distance_pixels")] = mean_distance_list
    mappings[("qc", "mean_centroid_height")] = mean_height_list
    mappings[("qc", "mean_centroid_width")] = mean_width_list

    # Spatial sorting makes the assigned IDs deterministic and makes nearby
    # cells appear near one another in the contour contact sheet.
    mappings = mappings.sort_values(
        by=[("qc", "mean_centroid_height"), ("qc", "mean_centroid_width")],
        kind="mergesort",
    ).copy()
    mappings.index = pd.Index(
        np.arange(1, len(mappings) + 1), name="global_cell_id"
    )

    validate_mapping_integrity(mappings, cents)
    n_wide = int(
        (
            (mappings[("qc", "n_sessions")] >= 2)
            & (mappings[("qc", "max_centroid_distance_pixels")] >= PARAM_DIST_PIXELS)
        ).sum()
    )
    if n_wide:
        logger.warning(
            "%d multi-session mappings have a final centroid span >= %.2f pixels. "
            "This can occur through transitive links; inspect them in the contour PDF.",
            n_wide,
            PARAM_DIST_PIXELS,
        )
    logger.info(
        "Mapping integrity check passed: every usable (session, unit_id) occurs once"
    )
    logger.info(
        "Assigned consecutive global_cell_id values 1..%d", len(mappings)
    )
    return mappings


def extract_matched_cell_contours(
    A_shifted,
    common_window_2d,
    mappings,
    cents,
    session_name_map,
    logger,
):
    """Extract compact vector contours for matched cells only."""
    session_cols = get_session_columns(mappings)
    cent_lookup = cents.set_index(["session", "unit_id"])[["height", "width"]]
    height_coords = np.asarray(A_shifted.coords["height"].values, dtype=float)
    width_coords = np.asarray(A_shifted.coords["width"].values, dtype=float)
    common_mask = np.asarray(
        common_window_2d.transpose("height", "width").compute().values,
        dtype=bool,
    )

    matched = mappings[mappings[("qc", "n_sessions")] >= CONTOUR_MIN_SESSIONS]
    logger.info(
        "Extracting contours for %d mappings present in at least %d sessions",
        len(matched),
        CONTOUR_MIN_SESSIONS,
    )
    records = []
    missing_contours = []
    processed_footprints = 0

    for global_cell_id, row in matched.iterrows():
        for col in session_cols:
            if pd.isna(row[col]):
                continue
            internal_session = col[1]
            display_session = session_name_map.get(internal_session, internal_session)
            unit_id = normalise_unit_id(row[col])
            footprint = (
                A_shifted.sel(session=internal_session, unit_id=unit_id)
                .transpose("height", "width")
                .compute()
            )
            values = np.asarray(footprint.values, dtype=float)
            values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
            values[~common_mask] = 0.0
            peak = float(values.max())
            if peak <= 0:
                missing_contours.append((global_cell_id, display_session, unit_id))
                continue

            level = peak * CONTOUR_LEVEL_FRACTION
            contours = [c for c in find_contours(values, level=level) if len(c) >= 3]
            if not contours:
                missing_contours.append((global_cell_id, display_session, unit_id))
                continue
            if CONTOUR_KEEP_LARGEST_ONLY:
                contours = [max(contours, key=len)]

            cent = cent_lookup.loc[(internal_session, unit_id)]
            centroid_y = float(cent["height"])
            centroid_x = float(cent["width"])
            for contour_id, contour in enumerate(contours):
                ys = np.interp(contour[:, 0], np.arange(len(height_coords)), height_coords)
                xs = np.interp(contour[:, 1], np.arange(len(width_coords)), width_coords)
                for point_order, (x_value, y_value) in enumerate(zip(xs, ys)):
                    records.append(
                        {
                            "global_cell_id": int(global_cell_id),
                            "session": display_session,
                            "unit_id": unit_id,
                            "contour_id": contour_id,
                            "point_order": point_order,
                            "x": round(float(x_value), 3),
                            "y": round(float(y_value), 3),
                            "centroid_x": round(centroid_x, 3),
                            "centroid_y": round(centroid_y, 3),
                            "contour_level": float(CONTOUR_LEVEL_FRACTION),
                        }
                    )
            processed_footprints += 1
            if processed_footprints % 100 == 0:
                logger.info("Extracted contours from %d footprints", processed_footprints)

    columns = [
        "global_cell_id",
        "session",
        "unit_id",
        "contour_id",
        "point_order",
        "x",
        "y",
        "centroid_x",
        "centroid_y",
        "contour_level",
    ]
    contour_df = pd.DataFrame.from_records(records, columns=columns)
    if missing_contours:
        logger.warning(
            "Could not extract %d contours; first examples: %s",
            len(missing_contours),
            missing_contours[:10],
        )
    return contour_df


def save_contour_contact_sheet(
    contour_df,
    mappings,
    session_name_map,
    figure_path,
    logger,
):
    """Save one small vector PDF with tiled local contour overlays."""
    if contour_df.empty:
        logger.warning("No matched-cell contours were available; no contour PDF saved")
        return

    display_sessions = list(session_name_map.values())
    cmap = plt.get_cmap("tab20")
    colors = {
        session: cmap(idx % cmap.N) for idx, session in enumerate(display_sessions)
    }
    handles = [
        Line2D([0], [0], color=colors[session], lw=1.5, label=session)
        for session in display_sessions
    ]
    # Show the widest (and therefore most suspicious) mappings first.
    global_ids = sorted(
        contour_df["global_cell_id"].unique(),
        key=lambda global_id: (
            -float(
                mappings.loc[
                    global_id, ("qc", "max_centroid_distance_pixels")
                ]
            ),
            int(global_id),
        ),
    )
    panels_per_page = CONTOUR_PANELS_PER_ROW * CONTOUR_ROWS_PER_PAGE
    figure_path = Path(figure_path)

    with PdfPages(figure_path) as pdf:
        for page_start in range(0, len(global_ids), panels_per_page):
            page_ids = global_ids[page_start : page_start + panels_per_page]
            fig, axes = plt.subplots(
                CONTOUR_ROWS_PER_PAGE,
                CONTOUR_PANELS_PER_ROW,
                figsize=(11.7, 8.3),
                squeeze=False,
            )
            axes_flat = axes.ravel()
            for ax, global_cell_id in zip(axes_flat, page_ids):
                cell_df = contour_df[
                    contour_df["global_cell_id"] == global_cell_id
                ]
                for (session, unit_id, contour_id), contour in cell_df.groupby(
                    ["session", "unit_id", "contour_id"], sort=False
                ):
                    ax.plot(
                        contour["x"],
                        contour["y"],
                        color=colors.get(session, "black"),
                        linewidth=0.9,
                    )
                    first = contour.iloc[0]
                    ax.plot(
                        first["centroid_x"],
                        first["centroid_y"],
                        marker="+",
                        color=colors.get(session, "black"),
                        markersize=4,
                        markeredgewidth=0.7,
                    )
                    ax.text(
                        first["centroid_x"],
                        first["centroid_y"],
                        str(int(unit_id)),
                        color=colors.get(session, "black"),
                        fontsize=4.5,
                        ha="left",
                        va="bottom",
                    )

                xmin, xmax = cell_df["x"].min(), cell_df["x"].max()
                ymin, ymax = cell_df["y"].min(), cell_df["y"].max()
                pad = CONTOUR_PADDING_PIXELS
                ax.set_xlim(xmin - pad, xmax + pad)
                ax.set_ylim(ymax + pad, ymin - pad)
                ax.set_aspect("equal", adjustable="box")
                ax.set_xticks([])
                ax.set_yticks([])
                max_dist = float(
                    mappings.loc[
                        global_cell_id, ("qc", "max_centroid_distance_pixels")
                    ]
                )
                n_sessions = int(
                    mappings.loc[global_cell_id, ("qc", "n_sessions")]
                )
                ax.set_title(
                    f"global {global_cell_id} | n={n_sessions} | max d={max_dist:.1f}px",
                    fontsize=6,
                )

            for ax in axes_flat[len(page_ids) :]:
                ax.axis("off")
            fig.suptitle(
                "Cross-registration QC: aligned footprint contours "
                f"({CONTOUR_LEVEL_FRACTION:.0%} of each footprint peak)",
                fontsize=10,
            )
            fig.legend(
                handles=handles,
                loc="upper center",
                bbox_to_anchor=(0.5, 0.965),
                ncol=min(5, max(1, len(handles))),
                fontsize=6,
                frameon=False,
            )
            fig.tight_layout(rect=(0.02, 0.02, 0.98, 0.91))
            pdf.savefig(fig)
            plt.close(fig)

    logger.info("Saved matched-cell contour contact sheet: %s", figure_path)


def restore_session_names(mappings, cents, shiftds, session_name_map):
    """Replace safe internal session labels with the original folder names."""
    mappings = mappings.copy()
    new_columns = []
    for col in mappings.columns:
        if isinstance(col, tuple) and col[0] == "session":
            new_columns.append((col[0], session_name_map.get(col[1], col[1])))
        else:
            new_columns.append(col)
    mappings.columns = pd.MultiIndex.from_tuples(new_columns)
    if ("group", "group") in mappings.columns:
        mappings = mappings.drop(columns=[("group", "group")])
    mappings = group_by_session(mappings)

    cents = cents.copy()
    cents["session"] = cents["session"].map(
        lambda session: session_name_map.get(session, session)
    )
    shiftds = shiftds.assign_coords(
        session=[
            session_name_map.get(session, session)
            for session in shiftds.coords["session"].values
        ]
    )
    return mappings, cents, shiftds


def count_registered_across_all_sessions(mappings_meta_fill):
    """Count mapping rows that contain a cell ID for every session."""
    session_cols = get_session_columns(mappings_meta_fill)

    registered_rows = mappings_meta_fill[
        mappings_meta_fill[session_cols].notna().all(axis=1)
    ]
    return int(len(registered_rows)), session_cols


def build_report_text(
    original_cells,
    original_cell_ids,
    cents,
    mappings,
    n_registered_across_all,
):
    """Build the final plain-text report that is appended to the log."""
    lines = []
    lines.append("")
    lines.append("=" * 80)
    lines.append("FINAL CROSS-REGISTRATION REPORT")
    lines.append("=" * 80)

    centroid_counts = cents.groupby("session")["unit_id"].nunique().to_dict()
    for session, n_cells in original_cells.items():
        unit_ids = np.asarray(original_cell_ids[session], dtype=float)
        if len(unit_ids):
            uid_min = int(unit_ids.min())
            uid_max = int(unit_ids.max())
            n_gaps = int(uid_max - uid_min + 1 - len(unit_ids))
            id_text = f"unit_id range = {uid_min}..{uid_max}; gaps = {n_gaps}"
        else:
            id_text = "no unit_id values"
        n_usable = int(centroid_counts.get(session, 0))
        lines.append(
            f"Session {session}: filtered C cells used = {n_cells}; {id_text}; "
            f"usable centroids in common FOV = {n_usable}; excluded = {n_cells - n_usable}"
        )

    lines.append("")
    lines.append(f"Global mapping rows = {len(mappings)}")
    lines.append(f"global_cell_id range = 1..{len(mappings)} (consecutive)")
    session_count_distribution = (
        mappings[("qc", "n_sessions")].value_counts().sort_index().to_dict()
    )
    for n_sessions, n_rows in session_count_distribution.items():
        lines.append(
            f"Rows containing cells from exactly {int(n_sessions)} session(s) = {int(n_rows)}"
        )
    lines.append(
        f"Cells registered across all sessions = {n_registered_across_all}"
    )
    lines.append(
        "Integrity check = passed (every usable session/unit_id appears exactly once)"
    )
    lines.append("=" * 80)
    return "\n".join(lines)


def build_mapping_csv_tables(mappings, c_column_lookups):
    """Build readable unit-ID and MATLAB C-column mapping tables."""
    session_cols = get_session_columns(mappings)
    unit_id_table = pd.DataFrame(
        {"global_cell_id": mappings.index.to_numpy(dtype=int)}
    )
    matlab_table = pd.DataFrame(
        {"global_cell_id": mappings.index.to_numpy(dtype=int)}
    )
    lookup_records = []

    for col in session_cols:
        session = col[1]
        if session not in c_column_lookups:
            raise ValueError(f"No C-column lookup found for session {session!r}")
        lookup = c_column_lookups[session]

        unit_ids = mappings[col].map(
            lambda value: pd.NA if pd.isna(value) else normalise_unit_id(value)
        ).astype("Int64")
        unit_id_table[session] = pd.array(unit_ids.to_numpy(), dtype="Int64")

        def to_matlab_column(value):
            if pd.isna(value):
                return pd.NA
            unit_id = normalise_unit_id(value)
            try:
                return lookup[unit_id]
            except KeyError as exc:
                raise ValueError(
                    f"Mapped unit_id {unit_id} from session {session!r} "
                    "has no corresponding column in C"
                ) from exc

        matlab_table[session] = pd.array(
            mappings[col].map(to_matlab_column).astype("Int64").to_numpy(),
            dtype="Int64",
        )
        lookup_records.extend(
            {
                "session": session,
                "unit_id": int(unit_id),
                "matlab_C_column": int(matlab_column),
            }
            for unit_id, matlab_column in lookup.items()
        )

    for qc_name in [
        "n_sessions",
        "max_centroid_distance_pixels",
        "mean_centroid_distance_pixels",
        "mean_centroid_height",
        "mean_centroid_width",
    ]:
        unit_id_table[qc_name] = mappings[("qc", qc_name)].to_numpy()

    lookup_table = pd.DataFrame.from_records(
        lookup_records,
        columns=["session", "unit_id", "matlab_C_column"],
    ).sort_values(["session", "matlab_C_column"], kind="mergesort")
    return unit_id_table, matlab_table, lookup_table


def save_outputs(
    dpath,
    mappings_meta_fill,
    cents,
    shiftds,
    contour_df,
    c_column_lookups,
    logger,
):
    """Save mapping, centroid, and shift outputs in DPATH."""
    dpath = Path(dpath)

    mappings_path = dpath / MAPPINGS_FILE_NAME
    mappings_csv_path = dpath / MAPPINGS_CSV_FILE_NAME
    matlab_mappings_csv_path = dpath / MATLAB_MAPPINGS_CSV_FILE_NAME
    unit_id_lookup_csv_path = dpath / UNIT_ID_LOOKUP_CSV_FILE_NAME
    cents_path = dpath / CENTROIDS_FILE_NAME
    shifts_path = dpath / SHIFTS_FILE_NAME
    contours_path = dpath / CONTOURS_FILE_NAME

    mappings_meta_fill.to_pickle(mappings_path)
    logger.info("Saved mappings: %s", mappings_path)

    unit_id_table, matlab_table, lookup_table = build_mapping_csv_tables(
        mappings_meta_fill,
        c_column_lookups,
    )
    unit_id_table.to_csv(mappings_csv_path, index=False)
    logger.info(
        "Saved CSV mappings with original MiniAn unit_id values: %s",
        mappings_csv_path,
    )
    matlab_table.to_csv(matlab_mappings_csv_path, index=False)
    logger.info(
        "Saved MATLAB-ready CSV mappings with 1-based C-column indices: %s",
        matlab_mappings_csv_path,
    )
    lookup_table.to_csv(unit_id_lookup_csv_path, index=False)
    logger.info("Saved unit_id to C-column lookup: %s", unit_id_lookup_csv_path)

    cents.to_pickle(cents_path)
    logger.info("Saved centroids: %s", cents_path)

    shiftds.to_netcdf(shifts_path)
    logger.info("Saved shifts dataset: %s", shifts_path)

    if SAVE_CONTOUR_QC:
        contour_df.to_csv(contours_path, index=False, compression="gzip")
        logger.info(
            "Saved compact contour coordinates: %s (%d coordinate points)",
            contours_path,
            len(contour_df),
        )


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Create one cross-session mapping for a single animal."
    )
    parser.add_argument(
        "--animal-root",
        required=True,
        help="Animal folder containing day/protocol/acquisition subfolders.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Optional output folder. By default a new timestamped folder is created under the animal folder.",
    )
    args = parser.parse_args()

    dpath = Path(args.animal_root).expanduser().resolve()
    if not dpath.is_dir():
        raise FileNotFoundError(
            "Animal root does not exist or is not a folder: {}".format(dpath)
        )

    if args.output_dir:
        output_dir = Path(args.output_dir).expanduser().resolve()
    else:
        output_dir = dpath / ("cross_registration_" + time.strftime("%Y%m%d_%H%M%S"))

    if output_dir == dpath:
        raise ValueError("Output directory must not be the animal root itself")

    log_file = output_dir / LOG_FILE_NAME
    logger = setup_logger(log_file)

    logger.info("Starting cross-session registration")
    logger.info("Animal root: %s", dpath)
    logger.info("Output directory: %s", output_dir)
    logger.info("Log file: %s", log_file)
    logger.info("PARAM_DIST_PIXELS: %s", PARAM_DIST_PIXELS)

    try:
        remove_previous_outputs(output_dir, logger)

        session_folders = find_session_miniscope_folders(dpath)
        if not session_folders:
            raise FileNotFoundError(
                f"No {MINISCOPE_FOLDER_NAME} folders found under {dpath}"
            )
        if len(session_folders) < 2:
            raise ValueError(
                "Cross-session registration needs at least two sessions. "
                f"Found only {len(session_folders)}."
            )

        unexpected_layout = [
            path
            for path in session_folders
            if len(path.relative_to(dpath).parts) != 3
        ]
        if unexpected_layout:
            raise ValueError(
                "Expected each MiniAn folder at "
                "ANIMAL_ROOT/day_or_protocol/acquisition/My_V4_Miniscope. "
                "Pass one animal folder as --animal-root, not the parent containing multiple animals. "
                f"First unexpected path: {unexpected_layout[0]}"
            )

        logger.info("Found %d sessions:", len(session_folders))
        for session_path in session_folders:
            logger.info("  %s", session_path)

        validate_minian_folders(session_folders)

        with ProgressBar(minimum=2):
            with timed_step(logger, "Opening MiniAn datasets"):
                (
                    minian_ds,
                    original_cells,
                    original_cell_ids,
                    session_name_map,
                    c_column_lookups,
                ) = open_all_sessions(session_folders, logger, dpath)

            with timed_step(logger, "Estimating session shifts"):
                shiftds = estimate_session_shifts(minian_ds)

            with timed_step(logger, "Creating common field-of-view window"):
                window, common_window_2d = make_common_window(shiftds)

            if HAS_MATPLOTLIB:
                save_alignment_figure(
                    shiftds,
                    common_window_2d,
                    output_dir / ALIGNMENT_FIGURE_NAME,
                    logger,
                    session_name_map=session_name_map,
                )
            else:
                logger.warning("matplotlib is unavailable; skipping alignment figure")

            A_shifted, cents, mappings_meta_fill = calculate_cross_session_mapping(
                minian_ds,
                shiftds,
                window,
                logger,
            )

            mappings_meta_fill = add_mapping_qc_and_global_ids(
                mappings_meta_fill,
                cents,
                logger,
            )

            if SAVE_CONTOUR_QC and HAS_MATPLOTLIB:
                with timed_step(logger, "Extracting compact matched-cell contours"):
                    contour_df = extract_matched_cell_contours(
                        A_shifted,
                        common_window_2d,
                        mappings_meta_fill,
                        cents,
                        session_name_map,
                        logger,
                    )
                save_contour_contact_sheet(
                    contour_df,
                    mappings_meta_fill,
                    session_name_map,
                    output_dir / CONTOUR_FIGURE_NAME,
                    logger,
                )
            else:
                contour_df = pd.DataFrame()
                if SAVE_CONTOUR_QC and not HAS_MATPLOTLIB:
                    logger.warning("matplotlib is unavailable; skipping contour QC figure")

        mappings_meta_fill, cents, shiftds = restore_session_names(
            mappings_meta_fill,
            cents,
            shiftds,
            session_name_map,
        )
        c_column_lookups = {
            session_name_map.get(session, session): lookup
            for session, lookup in c_column_lookups.items()
        }
        n_registered_across_all, session_cols = count_registered_across_all_sessions(
            mappings_meta_fill
        )
        logger.info("Mapping session columns: %s", session_cols)

        save_outputs(
            output_dir,
            mappings_meta_fill,
            cents,
            shiftds,
            contour_df,
            c_column_lookups,
            logger,
        )

        report_text = build_report_text(
            original_cells,
            original_cell_ids,
            cents,
            mappings_meta_fill,
            n_registered_across_all,
        )
        print(report_text)
        logger.info("\n%s", report_text)
        logger.info("Cross-session registration completed successfully")

    except Exception as exc:
        logger.error("Cross-session registration failed")
        logger.error(str(exc))
        logger.error(traceback.format_exc())
        raise


if __name__ == "__main__":
    main()
