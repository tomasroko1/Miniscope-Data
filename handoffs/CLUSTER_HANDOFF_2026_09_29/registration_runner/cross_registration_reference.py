"""Compare a new registration with the within-day mappings in merged MAT files."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io


def _session_key(path):
    parts = str(path).replace("\\", "/").strip("/").split("/")
    if len(parts) < 2:
        raise ValueError("Unexpected sess_paths entry: {!r}".format(path))
    return "/".join(parts[-2:])


def _mat_paths(session_struct):
    value = np.asarray(session_struct.sess_paths, dtype=object).reshape(-1)
    return [_session_key(item) for item in value]


def _matched_pairs(mapping, phase_a, phase_b):
    values = np.asarray(mapping[:, [phase_a, phase_b]], dtype=float)
    good = np.isfinite(values).all(axis=1)
    return {(int(a), int(b)) for a, b in values[good]}


def _new_pairs(mappings, session_a, session_b):
    columns = [("session", session_a), ("session", session_b)]
    values = mappings[columns].to_numpy(dtype=float)
    good = np.isfinite(values).all(axis=1)
    return {(int(a), int(b)) for a, b in values[good]}


def compare_to_mat_mappings(
    mappings, cents, reference_mat_dir, selected_sessions, output_dir, logger,
):
    """Compare exact MiniAn unit-ID pairs; global mapping row IDs are unrelated."""
    reference_mat_dir = Path(reference_mat_dir).expanduser().resolve()
    output_dir = Path(output_dir)
    if not reference_mat_dir.is_dir():
        raise FileNotFoundError("MAT reference folder not found: {}".format(reference_mat_dir))
    selected_sessions = set(selected_sessions)
    usable_by_session = {
        session: {int(value) for value in values}
        for session, values in cents.groupby("session")["unit_id"]
    }
    comparison_rows = []
    source_rows = []
    found_sessions = set()
    included_files = []

    for mat_path in sorted(reference_mat_dir.glob("*_merged.mat")):
        metadata = scipy.io.loadmat(
            str(mat_path), variable_names=["sess"], squeeze_me=True,
            struct_as_record=False,
        )
        if "sess" not in metadata or not hasattr(metadata["sess"], "sess_paths"):
            continue
        session_keys = _mat_paths(metadata["sess"])
        if not set(session_keys).issubset(selected_sessions):
            continue
        if len(session_keys) != 4:
            raise ValueError("{} does not have four reference acquisitions".format(mat_path))
        activity = scipy.io.loadmat(
            str(mat_path), variable_names=["act"], squeeze_me=True,
            struct_as_record=False,
        )
        if "act" not in activity or not hasattr(activity["act"], "mapping"):
            raise ValueError("{} has no act.mapping".format(mat_path))
        old_mapping = np.asarray(activity["act"].mapping, dtype=float)
        if old_mapping.ndim != 2 or old_mapping.shape[1] != len(session_keys):
            raise ValueError("Unexpected act.mapping shape in {}".format(mat_path))
        included_files.append(mat_path.name)
        found_sessions.update(session_keys)

        for phase, session in enumerate(session_keys):
            old_ids = {
                int(value) for value in old_mapping[:, phase]
                if np.isfinite(value)
            }
            new_ids = usable_by_session.get(session, set())
            shared = old_ids & new_ids
            source_rows.append({
                "reference_mat": mat_path.name,
                "session": session,
                "n_old_local_unit_ids": len(old_ids),
                "n_new_usable_unit_ids": len(new_ids),
                "n_unit_ids_in_both_sources": len(shared),
                "pct_old_ids_available_to_new": (
                    100.0 * len(shared) / len(old_ids) if old_ids else float("nan")
                ),
            })

        for phase_a in range(len(session_keys) - 1):
            for phase_b in range(phase_a + 1, len(session_keys)):
                session_a = session_keys[phase_a]
                session_b = session_keys[phase_b]
                old_pairs = _matched_pairs(old_mapping, phase_a, phase_b)
                new_pairs = _new_pairs(mappings, session_a, session_b)
                agreement = old_pairs & new_pairs
                n_old_a = int(np.isfinite(old_mapping[:, phase_a]).sum())
                n_old_b = int(np.isfinite(old_mapping[:, phase_b]).sum())
                n_new_a = len(usable_by_session.get(session_a, set()))
                n_new_b = len(usable_by_session.get(session_b, set()))
                comparison_rows.append({
                    "reference_mat": mat_path.name,
                    "session_a": session_a,
                    "session_b": session_b,
                    "n_old_matched_pairs": len(old_pairs),
                    "n_new_matched_pairs": len(new_pairs),
                    "n_identical_pairs": len(agreement),
                    "n_old_only_pairs": len(old_pairs - new_pairs),
                    "n_new_only_pairs": len(new_pairs - old_pairs),
                    "pct_old_pairs_recovered": (
                        100.0 * len(agreement) / len(old_pairs)
                        if old_pairs else float("nan")
                    ),
                    "pct_new_pairs_supported_by_old": (
                        100.0 * len(agreement) / len(new_pairs)
                        if new_pairs else float("nan")
                    ),
                    "n_old_cells_a": n_old_a,
                    "n_old_cells_b": n_old_b,
                    "n_new_usable_centroids_a": n_new_a,
                    "n_new_usable_centroids_b": n_new_b,
                    "old_match_pct_a": (
                        100.0 * len(old_pairs) / n_old_a if n_old_a else float("nan")
                    ),
                    "old_match_pct_b": (
                        100.0 * len(old_pairs) / n_old_b if n_old_b else float("nan")
                    ),
                    "new_match_pct_a": (
                        100.0 * len(new_pairs) / n_new_a if n_new_a else float("nan")
                    ),
                    "new_match_pct_b": (
                        100.0 * len(new_pairs) / n_new_b if n_new_b else float("nan")
                    ),
                })

    missing = sorted(selected_sessions - found_sessions)
    if missing:
        raise ValueError(
            "Selected acquisitions missing from reference MAT files: {}".format(missing)
        )
    if not comparison_rows:
        raise ValueError("No four-phase reference MAT files matched the selected sessions")
    pd.DataFrame(source_rows).to_csv(
        output_dir / "reference_unit_id_coverage.csv", index=False
    )
    pd.DataFrame(comparison_rows).to_csv(
        output_dir / "within_day_mapping_comparison.csv", index=False
    )
    report = {
        "reference_mat_files": included_files,
        "n_selected_sessions": len(selected_sessions),
        "n_within_day_pairs": len(comparison_rows),
        "note": (
            "Exact unit-ID pairs are compared. Within-day MAT mapping rows and "
            "new global_cell_id values are separate numbering systems. "
            "Percentages use their stated per-source denominators."
        ),
    }
    (output_dir / "within_day_mapping_comparison.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    logger.info(
        "Compared %d within-day acquisition pairs with MAT mappings: %s",
        len(comparison_rows), output_dir / "within_day_mapping_comparison.csv",
    )
    return comparison_rows

