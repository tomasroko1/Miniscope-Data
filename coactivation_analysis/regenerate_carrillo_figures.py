"""Regenerate Carrillo-Reid figures from saved matrices, without analysis."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from coactivation_analysis.carrillo_reid_pipeline import _plot_condition_clean


def regenerate(matrix_file: Path) -> list[Path]:
    output_dir = matrix_file.parent
    summary_file = output_dir / "summary.json"
    summary = json.loads(summary_file.read_text(encoding="utf-8")) if summary_file.exists() else {}
    with np.load(matrix_file, allow_pickle=False) as saved:
        arrays = {name: saved[name] for name in (
            "binary_activity", "similarity_binary", "singular_values",
            "state_labels", "significant_frames",
        )}

    phase_ranges = summary.get("plot_phase_ranges")
    phase_names = summary.get("plot_phase_names")
    phase_durations = summary.get("plot_phase_durations_seconds")
    condition = summary.get("condition", output_dir.name)
    figure = output_dir / "carrillo_reid_figure.png"
    _plot_condition_clean(
        arrays["binary_activity"], arrays["similarity_binary"],
        arrays["singular_values"], arrays["state_labels"], condition, figure,
        significant_frames=arrays["significant_frames"],
        phase_ranges=phase_ranges, phase_names=phase_names,
        phase_durations_seconds=phase_durations,
    )
    generated = [figure]
    vector_figure = output_dir / "carrillo_reid_figure_vector_order.png"
    if vector_figure.exists():
        if not (phase_ranges and phase_names):
            print(f"Skipped {vector_figure}: saved phase metadata is unavailable")
        else:
            _plot_condition_clean(
                arrays["binary_activity"], arrays["similarity_binary"],
                arrays["singular_values"], arrays["state_labels"], condition,
                vector_figure, significant_frames=arrays["significant_frames"],
                phase_ranges=phase_ranges, phase_names=phase_names,
                phase_durations_seconds=phase_durations, d_axis_mode="vector_order",
            )
            generated.append(vector_figure)
    return generated


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", type=Path, default=[Path("coactivation_analysis/results")])
    args = parser.parse_args()
    matrix_files = sorted({file for path in args.paths for file in (
        [path] if path.is_file() and path.name == "carrillo_reid_matrices.npz"
        else path.rglob("carrillo_reid_matrices.npz") if path.is_dir() else []
    )})
    if not matrix_files:
        parser.error("No carrillo_reid_matrices.npz files found in the supplied paths")
    for matrix_file in matrix_files:
        for figure in regenerate(matrix_file):
            print(figure)


if __name__ == "__main__":
    main()
