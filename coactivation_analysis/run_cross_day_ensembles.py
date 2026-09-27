"""Run Carrillo-Reid/Stoixeion on cells registered across different days.

The input ``mappings.csv`` is produced by ``run_cross_registration_by_animal.py``.
Each analysis compares a pair of days and a pair of phases, using only the
global cells that have a registered unit_id in both recordings.  The two
recordings are concatenated for one ensemble fit; the saved timeline marks the
boundary and labels both intervals.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path, PurePosixPath

import numpy as np
import scipy.io

if __package__:
    from .carrillo_reid_pipeline import (
        binarize_from_calcium_derivative,
        run_condition,
    )
else:
    from carrillo_reid_pipeline import (
        binarize_from_calcium_derivative,
        run_condition,
    )


STANDARD_PHASES = ("OF1", "SAMPLE", "TEST", "OF2")
HABL_PHASES = ("OF1", "OF2", "OF3", "OF4")


def _matlab_phases(value: object) -> list[np.ndarray]:
    if isinstance(value, np.ndarray) and value.dtype == object:
        return [np.asarray(item) for item in value.reshape(-1)]
    return [np.asarray(value)]


def _normalise_path(value: object) -> str:
    text = str(value).replace("\\", "/").strip("/")
    parts = [part for part in PurePosixPath(text).parts if part not in {"", "/"}]
    if len(parts) < 2:
        raise ValueError(f"No puedo identificar día/adquisición en {value!r}")
    return "/".join(parts[-2:])


def _scalar_text(value: object) -> str:
    if isinstance(value, np.ndarray):
        if value.size == 0:
            return ""
        if value.dtype == object:
            value = value.reshape(-1)[0]
        elif value.size == 1:
            value = value.item()
    if isinstance(value, bytes):
        value = value.decode(errors="replace")
    return str(value).strip()


def _mapping_lookup(mapping_csv: Path) -> tuple[dict[str, dict[int, int]], list[str]]:
    with mapping_csv.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "global_cell_id" not in reader.fieldnames:
            raise ValueError(f"{mapping_csv} no tiene columna global_cell_id")
        session_columns = [
            column for column in reader.fieldnames
            if column not in {
                "global_cell_id", "n_sessions", "max_centroid_distance_pixels",
                "mean_centroid_distance_pixels", "mean_centroid_height", "mean_centroid_width",
            }
        ]
        lookup: dict[str, dict[int, int]] = {}
        for row in reader:
            global_id = int(row["global_cell_id"])
            for column in session_columns:
                raw_unit_id = row.get(column, "")
                if raw_unit_id is None or not raw_unit_id.strip():
                    continue
                try:
                    unit_id = int(float(raw_unit_id))
                except ValueError as exc:
                    raise ValueError(
                        f"unit_id inválido en {column}: {raw_unit_id!r}"
                    ) from exc
                key = _normalise_path(column)
                if unit_id in lookup.setdefault(key, {}):
                    raise ValueError(f"unit_id repetido para {key}: {unit_id}")
                lookup[key][unit_id] = global_id
    return lookup, session_columns


def _load_day(mat_path: Path, animal: str, lookup: dict[str, dict[int, int]]) -> dict[str, object]:
    loaded = scipy.io.loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    if "act" not in loaded or "sess" not in loaded:
        raise ValueError("falta act o sess")
    act, sess = loaded["act"], loaded["sess"]
    if not hasattr(act, "S") or not hasattr(act, "C") or not hasattr(act, "t"):
        raise ValueError("falta act.S, act.C o act.t")
    if not hasattr(act, "mapping"):
        raise ValueError("act.mapping no está; no se pueden enlazar unit_id con cross-day")

    signals = {name: _matlab_phases(getattr(act, name)) for name in ("C", "S", "t")}
    mapping = np.asarray(act.mapping, dtype=float)
    if mapping.ndim == 1:
        mapping = mapping[:, None]
    paths_value = getattr(sess, "sess_paths", None)
    if paths_value is None:
        raise ValueError("sess.sess_paths no está")
    paths = [_scalar_text(item) for item in np.asarray(paths_value, dtype=object).reshape(-1)]
    n_phases = len(signals["S"])
    if len(paths) != n_phases or mapping.shape[1] != n_phases:
        raise ValueError(
            f"cantidad desigual de fases: S={n_phases}, paths={len(paths)}, mapping={mapping.shape[1]}"
        )

    day_name = _scalar_text(getattr(sess, "day_name", mat_path.stem))
    phase_names = HABL_PHASES if day_name.lower() == "habl" else STANDARD_PHASES
    if n_phases != len(phase_names):
        raise ValueError(f"se esperaban 4 fases para {day_name}; hay {n_phases}")

    records = []
    for phase_index in range(n_phases):
        session_key = _normalise_path(paths[phase_index])
        phase_lookup = lookup.get(session_key)
        if phase_lookup is None:
            raise KeyError(f"La adquisición {session_key!r} no está en mappings.csv")

        unit_ids = mapping[:, phase_index]
        valid_rows = np.flatnonzero(np.isfinite(unit_ids))
        ordered_rows = valid_rows[np.argsort(unit_ids[valid_rows])]
        ordered_ids = unit_ids[ordered_rows]
        C = np.asarray(signals["C"][phase_index], dtype=float)
        S = np.asarray(signals["S"][phase_index], dtype=float)
        t = np.asarray(signals["t"][phase_index], dtype=float).reshape(-1)
        if C.ndim == 1:
            C = C[:, None]
        if S.ndim == 1:
            S = S[:, None]
        if C.shape != S.shape or C.shape[0] != t.size or C.shape[1] != len(ordered_ids):
            raise ValueError(
                f"C/S/t/mapping no alinean en {session_key}: C={C.shape}, S={S.shape}, "
                f"t={t.shape}, IDs={len(ordered_ids)}"
            )
        if len(np.unique(ordered_ids)) != len(ordered_ids):
            raise ValueError(f"unit_id duplicados en act.mapping para {session_key}")

        by_global_id: dict[int, int] = {}
        missing = 0
        for column, raw_unit_id in enumerate(ordered_ids):
            unit_id = int(raw_unit_id)
            global_id = phase_lookup.get(unit_id)
            if global_id is None:
                missing += 1
            else:
                by_global_id[global_id] = column
        records.append({
            "animal": animal,
            "day_file": mat_path.name,
            "day_name": day_name,
            "phase": phase_names[phase_index],
            "session_key": session_key,
            "C": C,
            "S": S,
            "t": t,
            "by_global_id": by_global_id,
            "unmatched_unit_ids": missing,
        })
    return {"animal": animal, "day_file": mat_path.name, "day_name": day_name, "phases": records}


def _duration_seconds(t: np.ndarray) -> float:
    t = np.asarray(t, dtype=float)
    t = t[np.isfinite(t)]
    if t.size < 2:
        return 0.0
    steps = np.diff(t)
    steps = steps[np.isfinite(steps) & (steps > 0)]
    dt = float(np.median(steps)) if steps.size else 0.05
    return float(t.max() - t.min() + dt)


def _slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def run_pair(
    left: dict[str, object],
    right: dict[str, object],
    output_root: Path,
    *,
    minimum_cells: int,
    threshold_sd: float,
    pks_shuffles: int,
    scut_shuffles: int,
    hcut: float,
    seed: int,
) -> dict[str, object]:
    left_map = left["by_global_id"]
    right_map = right["by_global_id"]
    common_ids = np.asarray(sorted(set(left_map) & set(right_map)), dtype=int)
    label = f"{left['day_name']} {left['phase']} → {right['day_name']} {right['phase']}"
    result_row: dict[str, object] = {
        "animal": left["animal"],
        "source_day_file": left["day_file"],
        "source_phase": left["phase"],
        "source_acquisition": left["session_key"],
        "target_day_file": right["day_file"],
        "target_phase": right["phase"],
        "target_acquisition": right["session_key"],
        "n_shared_global_cells": int(len(common_ids)),
        "status": "skipped_too_few_shared_cells",
        "n_frames": 0,
        "n_significant_vectors": 0,
        "n_ensembles": 0,
        "pks": "",
        "scut": "",
        "output_dir": "",
    }
    if len(common_ids) < minimum_cells:
        return result_row

    left_columns = np.asarray([left_map[int(cell)] for cell in common_ids], dtype=int)
    right_columns = np.asarray([right_map[int(cell)] for cell in common_ids], dtype=int)
    left_events, left_thresholds, _ = binarize_from_calcium_derivative(
        left["C"][:, left_columns], left["S"][:, left_columns], left["t"], threshold_sd
    )
    right_events, right_thresholds, _ = binarize_from_calcium_derivative(
        right["C"][:, right_columns], right["S"][:, right_columns], right["t"], threshold_sd
    )
    events = np.vstack([left_events, right_events])
    thresholds = np.nanmean(np.vstack([left_thresholds, right_thresholds]), axis=0)
    ranges = [(0, len(left_events)), (len(left_events), len(events))]
    durations = [_duration_seconds(left["t"]), _duration_seconds(right["t"])]
    pair_dir = output_root / _slug(Path(str(left["day_file"])).stem) / (
        _slug(f"{left['phase']}_to_{Path(str(right['day_file'])).stem}_{right['phase']}")
    )
    fit = run_condition(
        None, None, None, common_ids, pair_dir,
        threshold_sd=threshold_sd,
        pks_shuffles=pks_shuffles,
        pks_null_vector_cap=1000,
        scut_shuffles=scut_shuffles,
        hcut=hcut,
        seed=seed,
        condition=label,
        events_override=events,
        event_thresholds_override=thresholds,
        event_definition=(
            f"separate per-recording dC/dt thresholds: positive dC/dt > "
            f"quiet-frame mean + {threshold_sd:g} SD; shared cells only"
        ),
        save_matrices=True,
        plot_phase_ranges=ranges,
        plot_phase_names=[
            f"{left['day_name']} · {left['phase']}",
            f"{right['day_name']} · {right['phase']}",
        ],
        plot_phase_durations_seconds=durations,
    )
    core_path = pair_dir / "core_neurons.csv"
    if core_path.exists():
        # run_condition uses a within-day column name; here IDs are from the
        # animal-wide cross-session registration.
        with core_path.open(newline="", encoding="utf-8") as handle:
            core_rows = list(csv.DictReader(handle))
        with core_path.open("w", newline="", encoding="utf-8") as handle:
            columns = ["condition", "factor", "rank", "global_cell_id"]
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            for core_row in core_rows:
                writer.writerow({
                    "condition": core_row.get("condition", label),
                    "factor": core_row["factor"],
                    "rank": core_row["rank"],
                    "global_cell_id": core_row["neuron_id_within_merged_day"],
                })
    result_row.update({
        "status": fit.get("status", "unknown"),
        "n_frames": fit.get("n_frames", len(events)),
        "n_significant_vectors": fit.get("n_significant_vectors", 0),
        "n_ensembles": fit.get("n_ensembles", 0),
        "pks": fit.get("pks", ""),
        "scut": fit.get("scut", ""),
        "output_dir": str(pair_dir),
    })
    (pair_dir / "cross_day_pair.json").write_text(
        json.dumps({**result_row, "shared_global_cell_ids": common_ids.tolist()}, indent=2),
        encoding="utf-8",
    )
    return result_row


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Ensemble analysis across days using cross-registered cells"
    )
    parser.add_argument("--animal", required=True, help="Animal label, e.g. R005")
    parser.add_argument("--data-dir", required=True, help="Folder containing DATA_DIR/ANIMAL/*.mat")
    parser.add_argument("--mappings-csv", required=True, help="mappings.csv from cross-registration")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--minimum-cells", type=int, default=20)
    parser.add_argument("--threshold-sd", type=float, default=3.0)
    parser.add_argument("--pks-shuffles", type=int, default=20)
    parser.add_argument("--scut-shuffles", type=int, default=20)
    parser.add_argument("--hcut", type=float, default=0.28)
    parser.add_argument("--seed", type=int, default=20260927)
    args = parser.parse_args()

    data_dir = Path(args.data_dir).expanduser().resolve() / args.animal
    output_dir = Path(args.output_dir).expanduser().resolve() / args.animal
    mapping_csv = Path(args.mappings_csv).expanduser().resolve()
    if not data_dir.is_dir():
        raise FileNotFoundError(f"No existe la carpeta de datos: {data_dir}")
    if not mapping_csv.is_file():
        raise FileNotFoundError(f"No existe mappings.csv: {mapping_csv}")
    output_dir.mkdir(parents=True, exist_ok=True)

    lookup, source_sessions = _mapping_lookup(mapping_csv)
    print(f"Mapping: {mapping_csv} ({len(source_sessions)} adquisiciones)")
    days = []
    skipped_days = []
    for mat_path in sorted(data_dir.glob("*.mat")):
        try:
            day = _load_day(mat_path, args.animal, lookup)
        except Exception as exc:
            skipped_days.append({"file": mat_path.name, "reason": str(exc)})
            print(f"[OMITIDA] {mat_path.name}: {exc}", flush=True)
            continue
        days.append(day)
        print(f"[DÍA] {mat_path.name}: {day['day_name']} ({len(day['phases'])} fases)", flush=True)

    # Compare matched ordinal subsessions across dates, including all four OF
    # sessions on HabL days. Add TEST -> OF2 only when those labels are real.
    phase_pairs = [(index, index) for index in range(4)]
    pair_rows = []
    for left_index, left_day in enumerate(days):
        for right_day in days[left_index + 1 :]:
            comparisons = list(phase_pairs)
            if (
                left_day["phases"][2]["phase"] == "TEST"
                and right_day["phases"][3]["phase"] == "OF2"
            ):
                comparisons.append((2, 3))
            for left_phase_index, right_phase_index in comparisons:
                left_phase = left_day["phases"][left_phase_index]
                right_phase = right_day["phases"][right_phase_index]
                row = run_pair(
                    left_phase,
                    right_phase,
                    output_dir,
                    minimum_cells=args.minimum_cells,
                    threshold_sd=args.threshold_sd,
                    pks_shuffles=args.pks_shuffles,
                    scut_shuffles=args.scut_shuffles,
                    hcut=args.hcut,
                    seed=args.seed + len(pair_rows),
                )
                pair_rows.append(row)
                print(
                    f"[{row['status']}] {row['source_day_file']} {row['source_phase']} → "
                    f"{row['target_day_file']} {row['target_phase']}: "
                    f"{row['n_shared_global_cells']} cells; {row['n_ensembles']} ensembles",
                    flush=True,
                )

    summary_path = output_dir / "cross_day_pairs.csv"
    if pair_rows:
        with summary_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(pair_rows[0]))
            writer.writeheader()
            writer.writerows(pair_rows)
    (output_dir / "cross_day_run.json").write_text(
        json.dumps({
            "animal": args.animal,
            "analysis": "pairwise cross-day Stoixeion on cross-registered cells",
            "mapping_csv": str(mapping_csv),
            "days_read": [day["day_file"] for day in days],
            "source_mapping_sessions": len(source_sessions),
            "minimum_cells": args.minimum_cells,
            "phase_pairs": ["same_subsession_index_across_days", "TEST:OF2_when_both_labels_exist"],
            "skipped_days": skipped_days,
            "n_pairwise_analyses": len(pair_rows),
            "n_analyses_run": sum(not str(row["status"]).startswith("skipped_") for row in pair_rows),
            "note": "Each fit uses only cells linked in both acquisitions; no cross-day cells are inferred from local column order.",
        }, indent=2),
        encoding="utf-8",
    )
    print(f"Resumen: {summary_path}")
    print(f"Resultados: {output_dir}")


if __name__ == "__main__":
    main()
