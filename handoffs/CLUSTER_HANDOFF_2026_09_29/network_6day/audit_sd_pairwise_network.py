"""Test same-day CA1 coactivity-network persistence in six SD sessions.

The network is a neuron-by-neuron correlation matrix of 1 s event-frame
fractions. This complementary analysis uses no object coordinates, position,
distance, or dwell. Animals, not neurons or edges, are treatment replicates.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse

from analyze_sd_object_code import SD_DAY_SPECS
from core.data_loader import get_common_neurons, get_data_dir, load_session, mapping_columns
from core.place_fields import threshold_event_frames


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "resultados_correlacion"
PHASES = ("OF1", "SAMPLE", "TEST", "OF2")
PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))


def _write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _plot_results(day_rows: list[dict], contrast_rows: list[dict], path: Path) -> None:
    day_rows = [
        row for row in day_rows
        if float(row["bin_width_s"]) == 1.0
        and row["global_mode"] == "population_mode_removed"
    ]
    contrast_rows = [
        row for row in contrast_rows
        if float(row["bin_width_s"]) == 1.0
        and row["global_mode"] == "population_mode_removed"
    ]
    if len(day_rows) != 6 or len(contrast_rows) != 3:
        raise ValueError("expected six day rows and three paired contrasts")
    colors = {"VEH": "#e07a1f", "CNO": "#3465a4"}
    animals = ("R004", "R005", "R006")
    fig, axes = plt.subplots(1, 2, figsize=(11.7, 4.7))
    ax = axes[0]
    for index, animal in enumerate(animals):
        for treatment, offset, marker in (("VEH", -0.12, "o"), ("CNO", 0.12, "s")):
            row = next(r for r in day_rows if r["animal"] == animal and r["treatment"] == treatment)
            x = index + offset
            observed = float(row["sample_test_full_observed"])
            null_q95 = float(row["sample_test_identity_null_q95"])
            ax.vlines(x, null_q95, observed, color=colors[treatment], alpha=0.65, lw=1.4)
            ax.scatter(x, observed, marker=marker, color=colors[treatment], s=58,
                       label=treatment if index == 0 else None, zorder=3)
            ax.scatter(x, null_q95, marker="_", color="#555555", s=180, zorder=3)
    ax.set_xticks(range(3), animals)
    ax.set_ylabel("Similitud de red SAMPLE–TEST")
    ax.set_title("Misma identidad neuronal vs. IDs permutados")
    ax.set_ylim(-0.03, 0.75)
    ax.grid(axis="y", alpha=0.2)
    ax.scatter([], [], marker="_", color="#555555", s=150, label="p95 IDs permutados")
    ax.legend(frameon=False, loc="upper left")

    ax = axes[1]
    for index, animal in enumerate(animals):
        row = next(r for r in contrast_rows if r["animal"] == animal)
        value = float(row["CNO_minus_VEH_sample_test_minus_of1_of2_mean"])
        ax.vlines(index, 0, value, color="#555555", lw=1.5)
        ax.scatter(index, value, color="#555555", s=65, zorder=3)
        ax.annotate(f"{value:+.2f}", (index, value), xytext=(0, 7 if value >= 0 else -14),
                    textcoords="offset points", ha="center", fontsize=9)
    ax.axhline(0, color="#777777", lw=1, ls="--")
    ax.set_xticks(range(3), animals)
    ax.set_ylabel("CNO−VEH del contraste\nSAMPLE–TEST menos OF1–OF2")
    ax.set_title("El contraste de tratamiento no es uniforme")
    ax.set_ylim(-0.42, 0.29)
    ax.grid(axis="y", alpha=0.2)
    fig.suptitle("Huella de coactividad de CA1 en seis jornadas S–D", y=0.99)
    id_p_max = max(float(row["sample_test_identity_null_p_upper"]) for row in day_rows)
    shift_p_max = max(float(row["sample_test_time_shift_null_p_upper"]) for row in day_rows)
    caption = (
        "Eventos S>3 SD, bins 1 s, modo poblacional removido; "
        + f"p de IDs≤{id_p_max:.3f} y p de corrimiento≤{shift_p_max:.3f} en 6/6 días. ".replace(".", ",")
        + "No es una prueba de CNO (n=3 animales)."
    )
    fig.text(0.5, 0.035, caption, ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.10, 1, 0.93), w_pad=3)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _binned_event_fractions(sub: dict, columns: np.ndarray, width_s: float,
                            event_mode: str = "frames") -> np.ndarray:
    t = np.asarray(sub["t"], dtype=float)
    s = np.asarray(sub["S"][:, columns], dtype=float)
    valid = np.isfinite(t) & np.isfinite(s).all(axis=1)
    t, s = t[valid], s[valid]
    if len(t) < 200 or not np.all(np.diff(t) > 0):
        raise ValueError("too few valid or non-monotone frames")
    event_matrix = np.column_stack([
        threshold_event_frames(s[:, col], threshold_sd=3.0)[0]
        for col in range(s.shape[1])
    ])
    if event_mode == "onsets":
        previous = np.vstack((np.zeros((1, event_matrix.shape[1])), event_matrix[:-1]))
        event_matrix = ((event_matrix > 0) & (previous == 0)).astype(float)
    elif event_mode != "frames":
        raise ValueError(event_mode)
    n_complete = int(np.floor((t[-1] - t[0]) / width_s))
    if n_complete < 40:
        raise ValueError("too few complete time bins")
    bin_index = np.floor((t - t[0]) / width_s).astype(int)
    keep = bin_index < n_complete
    bin_index, event_matrix = bin_index[keep], event_matrix[keep]
    incidence = sparse.csr_matrix(
        (np.ones(len(bin_index)), (bin_index, np.arange(len(bin_index)))),
        shape=(n_complete, len(bin_index)),
    )
    frames_per_bin = np.asarray(incidence @ np.ones(len(bin_index))).reshape(-1)
    typical = float(np.median(frames_per_bin[frames_per_bin > 0]))
    good = frames_per_bin >= 0.75 * typical
    if good.sum() < 40:
        raise ValueError("too few sufficiently sampled time bins")
    counts = np.asarray(incidence @ event_matrix)
    return counts[good] / frames_per_bin[good, None]


def _remove_population_mode(activity: np.ndarray) -> np.ndarray:
    """Regress each cell on the mean of all *other* cells in the same bins."""
    n_cells = activity.shape[1]
    if n_cells < 3:
        raise ValueError("population-mode control needs at least three cells")
    x = (activity.sum(axis=1, keepdims=True) - activity) / (n_cells - 1)
    y = activity
    xc = x - x.mean(axis=0, keepdims=True)
    yc = y - y.mean(axis=0, keepdims=True)
    denominator = np.sum(xc * xc, axis=0)
    beta = np.divide(
        np.sum(xc * yc, axis=0), denominator,
        out=np.zeros(n_cells), where=denominator > 0,
    )
    residual = yc - xc * beta[None, :]
    return residual


def _network_halves(day: dict, width_s: float, remove_population_mode: bool = False,
                    event_mode: str = "frames"):
    mapping = np.asarray(day["mapping"], dtype=float)
    common = get_common_neurons(mapping)
    if len(common) < 20:
        raise ValueError("too few identities mapped across four phases")
    binned = []
    for phase_index, sub in enumerate(day["subsessions"]):
        columns = mapping_columns(mapping, common, phase_index, sub["S"].shape[1])
        binned.append(_binned_event_fractions(sub, columns, width_s, event_mode))
    # Keep a cell only when it can contribute to both independent halves in
    # every phase. This prevents undefined pair correlations in one phase.
    eligibility = np.ones(len(common), dtype=bool)
    for phase in binned:
        for half in np.array_split(phase, 2):
            eligibility &= np.std(half, axis=0) > 0
            eligibility &= np.sum(half > 0, axis=0) >= 5
    common = common[eligibility]
    if len(common) < 20:
        raise ValueError("fewer than 20 informative cells after half-phase QC")
    matrices = []
    for phase in binned:
        half_matrices = []
        for half in np.array_split(phase, 2):
            activity = half[:, eligibility]
            if remove_population_mode:
                activity = _remove_population_mode(activity)
            half_matrices.append(np.corrcoef(activity, rowvar=False))
        matrices.append(tuple(half_matrices))
    if any(not np.isfinite(matrix).all() for pair in matrices for matrix in pair):
        raise AssertionError("non-finite pairwise correlation matrix")
    return common, matrices, [len(phase) for phase in binned], [
        phase[:, eligibility] for phase in binned
    ]


def _time_shift_test_null(matrices: list, binned: list[np.ndarray], width_s: float,
                          n_shuffles: int, rng, remove_population_mode: bool) -> np.ndarray:
    """Destroy synchronous TEST coactivity while preserving each cell's train."""
    test_halves = np.array_split(binned[2], 2)
    minimum = max(1, int(np.ceil(30.0 / width_s)))
    if any(len(half) <= 2 * minimum for half in test_halves):
        raise ValueError("TEST half too short for 30 s circular shifts")
    indices = np.arange(binned[2].shape[1])
    null = np.empty(n_shuffles, dtype=float)
    for shuffle_index in range(n_shuffles):
        shifted_matrices = []
        for half in test_halves:
            shifts = rng.integers(minimum, len(half) - minimum + 1, size=len(indices))
            shifted = np.column_stack([
                np.roll(half[:, col], int(shifts[col])) for col in indices
            ])
            if remove_population_mode:
                shifted = _remove_population_mode(shifted)
            shifted_matrices.append(np.corrcoef(shifted, rowvar=False))
        changed = list(matrices)
        changed[2] = tuple(shifted_matrices)
        null[shuffle_index] = _cross_phase(changed, 1, 2, indices)
    if not np.isfinite(null).all():
        raise AssertionError("non-finite time-shift null")
    return null


def _edge_similarity(a: np.ndarray, b: np.ndarray, indices: np.ndarray) -> float:
    i, j = np.triu_indices(len(indices), 1)
    aa = a[np.ix_(indices, indices)][i, j]
    bb = b[np.ix_(indices, indices)][i, j]
    if np.std(aa) == 0 or np.std(bb) == 0:
        return float("nan")
    return float(np.corrcoef(aa, bb)[0, 1])


def _cross_phase(matrices: list[tuple[np.ndarray, np.ndarray]], source: int,
                 target: int, indices: np.ndarray, permutation: np.ndarray | None = None) -> float:
    a0, a1 = matrices[source]
    b0, b1 = matrices[target]
    if permutation is not None:
        b0 = b0[np.ix_(permutation, permutation)]
        b1 = b1[np.ix_(permutation, permutation)]
    return float(np.mean((
        _edge_similarity(a0, b1, indices),
        _edge_similarity(a1, b0, indices),
    )))


def _day_metrics(matrices: list, indices: np.ndarray) -> dict[str, float]:
    metrics = {}
    for phase_index, phase_name in enumerate(PHASES):
        metrics[f"reliability_{phase_name}"] = _edge_similarity(
            matrices[phase_index][0], matrices[phase_index][1], indices
        )
    for source, target in PAIRS:
        metrics[f"similarity_{PHASES[source]}_{PHASES[target]}"] = _cross_phase(
            matrices, source, target, indices
        )
    metrics["sample_test_minus_of1_of2"] = (
        metrics["similarity_SAMPLE_TEST"] - metrics["similarity_OF1_OF2"]
    )
    return metrics


def _summarize(draws: list[dict[str, float]]) -> dict[str, float]:
    result = {}
    for key in draws[0]:
        values = np.asarray([row[key] for row in draws], dtype=float)
        result[f"{key}_mean"] = float(np.mean(values))
        result[f"{key}_q025"] = float(np.percentile(values, 2.5))
        result[f"{key}_q975"] = float(np.percentile(values, 97.5))
    return result


def _holm_adjusted(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    increasing = np.maximum.accumulate(
        (len(p_values) - np.arange(len(p_values))) * p_values[order]
    )
    adjusted = np.empty_like(p_values)
    adjusted[order] = np.minimum(1.0, increasing)
    return adjusted


def _self_test() -> None:
    rng = np.random.default_rng(612)
    n = 40
    latent = rng.normal(size=(500, 2))
    weights = rng.normal(size=(2, n))
    phase = latent @ weights + rng.normal(scale=0.4, size=(500, n))
    matrices = [tuple(np.corrcoef(half, rowvar=False) for half in np.array_split(phase, 2))] * 4
    indices = np.arange(n)
    matched = _cross_phase(matrices, 0, 1, indices)
    permuted = _cross_phase(matrices, 0, 1, indices, rng.permutation(n))
    assert matched > 0.9 and permuted < 0.3
    assert abs(_day_metrics(matrices, indices)["sample_test_minus_of1_of2"]) < 1e-12
    tiny = {
        "t": np.arange(1000) * 0.05,
        "S": rng.binomial(1, 0.05, size=(1000, 5)),
    }
    assert 48 <= len(_binned_event_fractions(tiny, np.arange(5), 1.0)) <= 50
    assert 48 <= len(_binned_event_fractions(tiny, np.arange(5), 1.0, "onsets")) <= 50
    pop = np.linspace(-1, 1, 500)[:, None] * rng.uniform(0.5, 2, size=(1, n))
    residual = _remove_population_mode(pop)
    assert np.max(np.abs(residual)) < 1e-10
    synthetic_binned = [phase] * 4
    shifted_null = _time_shift_test_null(matrices, synthetic_binned, 1.0, 5, rng, False)
    assert np.max(shifted_null) < matched - 0.4
    assert np.allclose(_holm_adjusted(np.repeat(0.005, 6)), 0.03)
    print("Network matrix, identity permutation, and time-binning checks passed.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path(get_data_dir()))
    parser.add_argument("--cell-subsamples", type=int, default=200)
    parser.add_argument("--identity-shuffles", type=int, default=199)
    parser.add_argument("--time-shuffles", type=int, default=199)
    parser.add_argument("--event-mode", choices=("frames", "onsets"), default="frames")
    parser.add_argument("--seed", type=int, default=6122026)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        _self_test()
        return
    if args.plot_only:
        suffix = "" if args.event_mode == "frames" else "_onsets"
        with (OUT / f"SD_pairwise_network_by_day{suffix}.csv").open(newline="", encoding="utf-8") as handle:
            day_rows = list(csv.DictReader(handle))
        with (OUT / f"SD_pairwise_network_contrasts{suffix}.csv").open(newline="", encoding="utf-8") as handle:
            contrast_rows = list(csv.DictReader(handle))
        if args.event_mode == "frames":
            _plot_results(day_rows, contrast_rows, OUT / "17_SD_pairwise_network.png")
        return
    days = [
        (animal, file, treatment, load_session(animal, file, data_dir=str(args.data_dir)))
        for animal, file, treatment in SD_DAY_SPECS
    ]
    rows = []
    contrasts = []
    for width_s in (1.0, 0.5, 2.0):
      for global_mode in ("raw", "population_mode_removed"):
        cache = {}
        for animal, file, treatment, day in days:
            ids, matrices, bin_counts, binned = _network_halves(
                day, width_s, remove_population_mode=global_mode != "raw",
                event_mode=args.event_mode,
            )
            cache[(animal, treatment)] = (ids, matrices, bin_counts, binned)
        for animal_index, animal in enumerate(("R004", "R005", "R006")):
            minimum_cells = min(len(cache[(animal, treatment)][0]) for treatment in ("VEH", "CNO"))
            draws_by_treatment = {}
            for treatment_index, treatment in enumerate(("VEH", "CNO")):
                ids, matrices, bin_counts, binned = cache[(animal, treatment)]
                width_key = int(round(100 * width_s))
                mode_key = 0 if global_mode == "raw" else 1
                subset_rng = np.random.default_rng(np.random.SeedSequence(
                    [args.seed, animal_index, treatment_index, width_key, 0]
                ))
                identity_rng = np.random.default_rng(np.random.SeedSequence(
                    [args.seed, animal_index, treatment_index, width_key, mode_key, 1]
                ))
                time_rng = np.random.default_rng(np.random.SeedSequence(
                    [args.seed, animal_index, treatment_index, width_key, mode_key, 2]
                ))
                subsets = [
                    np.sort(subset_rng.choice(len(ids), size=minimum_cells, replace=False))
                    for _ in range(args.cell_subsamples)
                ]
                draws = [_day_metrics(matrices, selected) for selected in subsets]
                draws_by_treatment[treatment] = draws
                full_indices = np.arange(len(ids))
                observed = _cross_phase(matrices, 1, 2, full_indices)
                identity_null = np.asarray([
                    _cross_phase(matrices, 1, 2, full_indices, identity_rng.permutation(len(ids)))
                    for _ in range(args.identity_shuffles)
                ])
                if not np.isfinite(identity_null).all():
                    raise AssertionError("non-finite identity null")
                if width_s == 1.0 and global_mode == "population_mode_removed":
                    time_null = _time_shift_test_null(
                        matrices, binned, width_s, args.time_shuffles, time_rng, True
                    )
                    time_null_q95 = float(np.percentile(time_null, 95))
                    time_null_p = float((1 + (time_null >= observed).sum()) / (len(time_null) + 1))
                else:
                    time_null_q95 = float("nan")
                    time_null_p = float("nan")
                rows.append({
                    "animal": animal,
                    "treatment": treatment,
                    "day_id": f"{animal}_{Path(next(file for a, file, tr, _ in days if a == animal and tr == treatment)).stem}",
                    "bin_width_s": width_s,
                    "event_mode": args.event_mode,
                    "global_mode": global_mode,
                    "n_cells_after_half_qc": len(ids),
                    "n_cells_equalized": minimum_cells,
                    "n_bins_OF1": bin_counts[0],
                    "n_bins_SAMPLE": bin_counts[1],
                    "n_bins_TEST": bin_counts[2],
                    "n_bins_OF2": bin_counts[3],
                    "sample_test_full_observed": observed,
                    "sample_test_identity_null_mean": float(np.mean(identity_null)),
                    "sample_test_identity_null_q95": float(np.percentile(identity_null, 95)),
                    "sample_test_identity_null_p_upper": float((1 + (identity_null >= observed).sum()) / (len(identity_null) + 1)),
                    "sample_test_time_shift_null_q95": time_null_q95,
                    "sample_test_time_shift_null_p_upper": time_null_p,
                    **_summarize(draws),
                })
            delta = np.asarray([
                c["sample_test_minus_of1_of2"] - v["sample_test_minus_of1_of2"]
                for v, c in zip(draws_by_treatment["VEH"], draws_by_treatment["CNO"])
            ])
            contrasts.append({
                "animal": animal,
                "bin_width_s": width_s,
                "event_mode": args.event_mode,
                "global_mode": global_mode,
                "n_cells_equalized": minimum_cells,
                "CNO_minus_VEH_sample_test_minus_of1_of2_mean": float(np.mean(delta)),
                "cell_subset_q025": float(np.percentile(delta, 2.5)),
                "cell_subset_q975": float(np.percentile(delta, 97.5)),
            })
            print(animal, width_s, global_mode,
                  f"CNO-VEH task-vs-return={np.mean(delta):+.3f}", flush=True)
    primary = [
        row for row in rows
        if row["bin_width_s"] == 1.0
        and row["global_mode"] == "population_mode_removed"
    ]
    if len(primary) != 6:
        raise AssertionError("Holm family must contain exactly six SD days")
    for source, target in (
        ("sample_test_identity_null_p_upper", "sample_test_identity_holm_p"),
        ("sample_test_time_shift_null_p_upper", "sample_test_time_shift_holm_p"),
    ):
        adjusted = _holm_adjusted(np.asarray([row[source] for row in primary]))
        for row, value in zip(primary, adjusted):
            row[target] = float(value)
    for row in rows:
        row.setdefault("sample_test_identity_holm_p", float("nan"))
        row.setdefault("sample_test_time_shift_holm_p", float("nan"))
    OUT.mkdir(parents=True, exist_ok=True)
    suffix = "" if args.event_mode == "frames" else "_onsets"
    _write_csv(OUT / f"SD_pairwise_network_by_day{suffix}.csv", rows)
    _write_csv(OUT / f"SD_pairwise_network_contrasts{suffix}.csv", contrasts)
    (OUT / f"SD_pairwise_network_protocol{suffix}.json").write_text(json.dumps({
        "method": "Pearson correlation of neuron-neuron edge matrices from within-phase binned event-frame fractions",
        "primary_bin_width_s": 1.0,
        "sensitivity_bin_width_s": [0.5, 2.0],
        "event_threshold": "S > 3*within-cell/phase SD(S)",
        "event_mode": args.event_mode,
        "event_unit": "threshold-positive frames" if args.event_mode == "frames" else "first frame of each threshold-positive run",
        "population_mode_sensitivity": "regress each cell on the mean of all other cells within the same half-phase before pairwise correlations",
        "cell_registration": "mapping common to four phases within each day only",
        "cell_qc": "nonconstant and at least five active bins in every phase half",
        "within_phase_split": "two time halves, independently estimated edge matrices",
        "cross_phase_similarity": "Pearson across upper-triangle edges, average of half1-to-half2 and half2-to-half1",
        "day_contrast": "SAMPLE-to-TEST network similarity minus OF1-to-OF2 network similarity",
        "treatment_contrast": "paired CNO minus VEH within animal, no between-date cell tracking",
        "identity_null": f"{args.identity_shuffles} random neuron-label permutations in TEST for SAMPLE-to-TEST matching; one-sided p unadjusted",
        "time_shift_null": f"for primary 1 s binned population-mode-removed test only, {args.time_shuffles} independent circular TEST shifts per cell/half by at least 30 s",
        "multiple_testing": "Holm adjustment separately across the six primary SD days for each of the identity and time-shift null families",
        "cell_subsamples": args.cell_subsamples,
        "seed": args.seed,
        "random_streams": "independent SeedSequence streams for cell subsets, identity null, and time-shift null; cell subsets shared across global-mode sensitivity",
        "inference_unit": "animal (n=3); cell-subset percentiles are not biological confidence intervals",
        "no_position_or_object_behavior": True,
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    if args.event_mode == "frames":
        _plot_results(rows, contrasts, OUT / "17_SD_pairwise_network.png")


if __name__ == "__main__":
    main()
