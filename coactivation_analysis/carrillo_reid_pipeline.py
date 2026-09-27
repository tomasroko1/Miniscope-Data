"""Python port of Carrillo-Reid's Stoixeion ensemble-identification steps.

The binary input is built from the first derivative of the continuous calcium
trace C.  Per-cell noise is estimated from frames where the supplied
deconvolution S is zero, and activity is set at baseline mean + 3 SD.  The
remaining stages follow Stoixeion.m: shuffle-based pks (optionally with the
source MATLAB histogram cutoff), TF-IDF, cosine thresholding, two Jaccard
denoising passes, and the SVD
factor reconstruction/cutoff used by Edos_from_Sindex_svd.m.

The MATLAB method uses a single recording-wide binary raster.  This module
accepts one subsession at a time so the four within-day conditions can be
compared using the cell registration that is actually present in each merged
file.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from scipy.stats import rankdata


@dataclass
class BinarizationInfo:
    threshold_sd: float
    active_fraction: float
    median_threshold: float
    median_noise_mean: float
    median_noise_sd: float
    cells_using_fallback: int


def _finite_sd(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.std(values, ddof=1)) if values.size > 1 else 0.0


def binarize_from_calcium_derivative(
    fluorescence: np.ndarray,
    deconvolved: np.ndarray,
    timestamps: np.ndarray,
    threshold_sd: float = 3.0,
) -> tuple[np.ndarray, np.ndarray, BinarizationInfo]:
    """Threshold positive fluorescence changes against quiescent-frame noise.

    ``S == 0`` marks candidate quiet samples.  The threshold is computed
    independently for each cell as mean(dC/dt | quiet) + threshold_sd times
    SD(dC/dt | quiet).  If a cell has fewer than 100 usable quiet samples, all
    finite derivatives are used and that fallback is reported.
    """
    c = np.asarray(fluorescence, dtype=float)
    s = np.asarray(deconvolved, dtype=float)
    t = np.asarray(timestamps, dtype=float)
    if c.ndim != 2 or s.shape != c.shape or len(t) != len(c):
        raise ValueError("C, S y t deben alinearse como frames x neuronas")
    if len(t) < 2:
        raise ValueError("Se necesitan al menos dos frames para calcular dC/dt")

    dt = np.diff(t)
    valid_dt = np.isfinite(dt) & (dt > 0)
    if not np.any(valid_dt):
        raise ValueError("No hay intervalos de tiempo positivos")
    typical_dt = float(np.median(dt[valid_dt]))
    dt = np.where(valid_dt, dt, typical_dt)
    derivative = np.diff(c, axis=0) / dt[:, None]
    event_raster = np.zeros(c.shape, dtype=np.uint8)
    thresholds = np.full(c.shape[1], np.nan, dtype=float)
    noise_means = np.full(c.shape[1], np.nan, dtype=float)
    noise_sds = np.full(c.shape[1], np.nan, dtype=float)
    fallbacks = 0

    for cell in range(c.shape[1]):
        d = derivative[:, cell]
        quiet = (s[1:, cell] <= 0) & np.isfinite(d)
        noise = d[quiet]
        if noise.size < 100:
            noise = d[np.isfinite(d)]
            fallbacks += 1
        if noise.size < 2:
            continue
        noise_mean = float(np.mean(noise))
        noise_sd = _finite_sd(noise)
        threshold = noise_mean + threshold_sd * noise_sd
        thresholds[cell] = threshold
        noise_means[cell] = noise_mean
        noise_sds[cell] = noise_sd
        event_raster[1:, cell] = (
            np.isfinite(d) & (d > threshold)
        ).astype(np.uint8)

    info = BinarizationInfo(
        threshold_sd=float(threshold_sd),
        active_fraction=float(event_raster.mean()),
        median_threshold=float(np.nanmedian(thresholds)),
        median_noise_mean=float(np.nanmedian(noise_means)),
        median_noise_sd=float(np.nanmedian(noise_sds)),
        cells_using_fallback=int(fallbacks),
    )
    return event_raster, thresholds, info


def _column_cosine_means(vectors: np.ndarray) -> np.ndarray:
    """Mean cosine of each column to all selected columns, diagonal included."""
    x = np.asarray(vectors, dtype=float)
    if x.ndim != 2 or x.shape[1] == 0:
        return np.empty(0, dtype=float)
    norms = np.linalg.norm(x, axis=0)
    valid = norms > 0
    if not np.any(valid):
        return np.empty(0, dtype=float)
    y = x[:, valid] / norms[valid][None, :]
    return (y.T @ y.sum(axis=1)) / y.shape[1]


def _time_shuffle_each_cell(
    activity_cells_by_frame: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    shuffled = np.empty_like(activity_cells_by_frame)
    for cell in range(activity_cells_by_frame.shape[0]):
        shuffled[cell] = activity_cells_by_frame[cell, rng.permutation(activity_cells_by_frame.shape[1])]
    return shuffled


def stoixeion_histc_cutoff(values: np.ndarray, percentile: float = 98.0) -> tuple[float, int]:
    """Reproduce the source MATLAB `histc(..., 0:0.02:max(...))` cutoff.

    The last listed bin counts only values equal to that edge; values above
    it are omitted by MATLAB's `histc`. This is a fidelity option, not a
    more rigorous null than the continuous empirical percentile.
    """
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    if not data.size:
        return float("nan"), 0
    edges = np.arange(0.0, float(np.max(data)) + 1e-12, 0.02)
    if not len(edges):
        return float("nan"), len(data)
    bins = np.searchsorted(edges, data, side="right") - 1
    included = (bins >= 0) & (bins < len(edges) - 1)
    included |= np.isclose(data, edges[-1], atol=1e-12, rtol=0.0)
    counts = np.bincount(bins[included], minlength=len(edges))
    if counts.sum() == 0:
        return float("nan"), len(data)
    cumulative = np.cumsum(counts) / counts.sum()
    selected = np.flatnonzero(cumulative > percentile / 100.0)
    cutoff = float(edges[selected[0]]) if len(selected) else float("nan")
    return cutoff, int((~included).sum())


def select_significant_vectors(
    events_frames_by_cells: np.ndarray,
    n_shuffles: int = 100,
    null_percentile: float = 98.0,
    seed: int = 1729,
    min_coactive_cells: int = 3,
    max_null_vectors_per_shuffle: int | None = 2000,
    cutoff_method: str = "continuous_percentile",
) -> tuple[np.ndarray, int, dict[str, float | int | str]]:
    """Port of Stoixeion.findHighactFrames with its automatic pks criterion."""
    if cutoff_method not in ("continuous_percentile", "stoixeion_histc"):
        raise ValueError(f"Unknown pks cutoff method: {cutoff_method}")
    events = np.asarray(events_frames_by_cells, dtype=np.uint8).T
    frame_counts = events.sum(axis=0)
    maximum = int(frame_counts.max(initial=0))
    if maximum < min_coactive_cells:
        return np.empty(0, dtype=int), maximum + 1, {
            "status": "no_frame_reached_minimum_coactivity",
            "null_percentile": float(null_percentile),
            "n_shuffles": int(n_shuffles),
            "cutoff_method": cutoff_method,
        }

    rng = np.random.default_rng(seed)
    shuffled = [
        _time_shuffle_each_cell(events, rng)
        for _ in range(n_shuffles)
    ]
    shuffled_counts = [x.sum(axis=0) for x in shuffled]

    diagnostic: dict[str, float | int | str] = {
        "status": "ok",
        "null_percentile": float(null_percentile),
        "n_shuffles": int(n_shuffles),
        "cutoff_method": cutoff_method,
        "max_null_vectors_per_shuffle": (
            "all"
            if max_null_vectors_per_shuffle is None or max_null_vectors_per_shuffle <= 0
            else int(max_null_vectors_per_shuffle)
        ),
    }
    selected_pks = maximum + 1
    selected_real_mean = np.nan
    selected_null_cutoff = np.nan
    selected_count = 0

    for pks in range(min_coactive_cells, maximum + 1):
        real_vectors = events[:, frame_counts >= pks]
        real_scores = _column_cosine_means(real_vectors)
        if real_scores.size < 2:
            continue
        null_scores = []
        for sh, counts in zip(shuffled, shuffled_counts):
            null_indices = np.flatnonzero(counts >= pks)
            if (
                max_null_vectors_per_shuffle is not None
                and max_null_vectors_per_shuffle > 0
                and null_indices.size > max_null_vectors_per_shuffle
            ):
                null_indices = np.sort(
                    rng.choice(null_indices, size=max_null_vectors_per_shuffle, replace=False)
                )
            null_vectors = sh[:, null_indices]
            scores = _column_cosine_means(null_vectors)
            if scores.size:
                null_scores.append(scores)
        if not null_scores:
            continue
        null_values = np.concatenate(null_scores)
        if cutoff_method == "stoixeion_histc":
            cutoff, excluded = stoixeion_histc_cutoff(null_values, null_percentile)
        else:
            cutoff = float(np.percentile(null_values, null_percentile))
            excluded = 0
        real_mean = float(np.mean(real_scores))
        if real_mean > cutoff:
            selected_pks = pks
            selected_real_mean = real_mean
            selected_null_cutoff = cutoff
            selected_count = int(real_vectors.shape[1])
            diagnostic["null_values_outside_histc_bins"] = int(excluded)
            break

    if selected_count == 0:
        diagnostic["status"] = "no_real_similarity_exceeded_shuffle_cutoff"
        return np.empty(0, dtype=int), selected_pks, diagnostic

    diagnostic.update(
        {
            "real_mean_similarity": selected_real_mean,
            "shuffle_similarity_cutoff": selected_null_cutoff,
            "selected_vectors": selected_count,
        }
    )
    return np.flatnonzero(frame_counts >= selected_pks), selected_pks, diagnostic


def select_vectors_at_fixed_pks(
    events_frames_by_cells: np.ndarray, pks: int, source: str
) -> tuple[np.ndarray, dict[str, float | int | str]]:
    """Apply the control-condition pks cutoff unchanged to a comparison segment."""
    counts = np.asarray(events_frames_by_cells, dtype=np.uint8).sum(axis=1)
    frames = np.flatnonzero(counts >= int(pks))
    return frames, {
        "status": "fixed_from_control_condition",
        "pks_source": source,
        "selected_vectors": int(len(frames)),
    }


def stoixeion_tfidf(population_vectors: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """TF-IDF weighting in the original Stoixeion Ras_tf_idf.m implementation."""
    binary = np.asarray(population_vectors, dtype=float)
    n_cells, n_vectors = binary.shape
    frame_activity = binary.sum(axis=0)
    tf = binary / np.where(frame_activity > 0, frame_activity, 1.0)[None, :]
    document_frequency = (binary == 1).sum(axis=1)
    idf = np.ones(n_cells, dtype=float)
    present = document_frequency > 0
    idf[present] += np.log(n_vectors / document_frequency[present])
    idf[~present] += np.log(max(n_vectors, 1))
    return tf * idf[:, None], idf


def cosine_similarity_map(vectors: np.ndarray) -> np.ndarray:
    x = np.asarray(vectors, dtype=float)
    norms = np.linalg.norm(x, axis=0, keepdims=True)
    normalized = x / np.where(norms > 0, norms, 1.0)
    similarity = np.clip(normalized.T @ normalized, 0.0, 1.0)
    np.fill_diagonal(similarity, 1.0)
    return similarity


def estimate_scut_from_shuffles(
    tfidf_vectors: np.ndarray,
    n_shuffles: int = 20,
    percentile: float = 88.0,
    seed: int = 2718,
    pair_samples: int = 100_000,
) -> float:
    """Stoixeion.calc_scut threshold; large maps use Monte Carlo pair sampling."""
    x = np.asarray(tfidf_vectors, dtype=float)
    n_states = x.shape[1]
    rng = np.random.default_rng(seed)
    bins = np.linspace(0.0, 1.0, 101)
    histogram = np.zeros(len(bins) - 1, dtype=np.int64)
    for _ in range(n_shuffles):
        shuffled = np.empty_like(x)
        for neuron in range(x.shape[0]):
            shuffled[neuron] = x[neuron, rng.permutation(n_states)]
        if n_states <= 1500:
            similarity = cosine_similarity_map(shuffled)
            counts, _ = np.histogram(similarity.ravel(), bins=bins)
        else:
            norms = np.linalg.norm(shuffled, axis=0)
            normalized = shuffled / np.where(norms > 0, norms, 1.0)[None, :]
            left = rng.integers(0, n_states, size=pair_samples)
            right = rng.integers(0, n_states, size=pair_samples)
            keep = left != right
            left, right = left[keep], right[keep]
            sampled = np.empty(len(left), dtype=np.float32)
            chunk = 10_000
            for start in range(0, len(left), chunk):
                stop = min(start + chunk, len(left))
                sampled[start:stop] = np.einsum(
                    "ij,ij->j",
                    normalized[:, left[start:stop]],
                    normalized[:, right[start:stop]],
                    optimize=True,
                )
            counts, _ = np.histogram(np.clip(sampled, 0.0, 1.0), bins=bins)
        histogram += counts
    cumulative = np.cumsum(histogram, dtype=float) / max(int(histogram.sum()), 1)
    index = int(np.flatnonzero(cumulative > percentile / 100.0)[0])
    return float(bins[index])


def _jaccard_similarity_between_rows(binary_rows: np.ndarray) -> np.ndarray:
    binary_rows = np.asarray(binary_rows, dtype=np.uint8)
    if binary_rows.shape[0] > 1500 and float(binary_rows.mean()) < 0.5:
        from scipy.sparse import csr_matrix

        sparse_rows = csr_matrix(binary_rows.astype(np.int32, copy=False))
        intersections = (sparse_rows @ sparse_rows.T).toarray()
        row_counts = np.asarray(sparse_rows.sum(axis=1)).ravel().astype(np.int64)
    else:
        x = binary_rows.astype(np.int32, copy=False)
        intersections = x @ x.T
        row_counts = x.sum(axis=1, dtype=np.int64)
    unions = row_counts[:, None] + row_counts[None, :] - intersections
    similarity = np.divide(
        intersections,
        unions,
        out=np.ones(intersections.shape, dtype=np.float32),
        where=unions > 0,
    )
    np.fill_diagonal(similarity, 1.0)
    return similarity


def two_pass_hamming_denoise(
    binary_similarity: np.ndarray, hcut: float
) -> np.ndarray:
    """Two Jaccard-similarity passes matching Stoixeion's `1-Hdist`.

    Despite its name, the MATLAB Hdist divides XOR by OR, so it is Jaccard
    distance rather than ordinary Hamming distance.
    """
    first_similarity = _jaccard_similarity_between_rows(binary_similarity)
    first_binary = (first_similarity > hcut).astype(np.uint8)
    del first_similarity
    second_similarity = _jaccard_similarity_between_rows(first_binary)
    final_binary = (second_similarity > hcut).astype(np.uint8)
    del second_similarity
    final_binary = np.maximum(final_binary, final_binary.T)
    np.fill_diagonal(final_binary, 1)
    return final_binary


def svd_ensemble_support(
    final_binary_map: np.ndarray,
    n_neurons: int,
    minimum_prevalence: float = 0.05,
    initial_factor_cut: float = 0.4,
    factor_cut_step: float = 0.01,
    return_singular_ranks: bool = False,
) -> tuple:
    """Port of Edos_from_Sindex_svd.m, including factor-cut iteration."""
    matrix = np.asarray(final_binary_map, dtype=np.float32)
    if matrix.ndim != 2 or matrix.shape[0] < 2:
        raise ValueError("El mapa final debe ser cuadrado y tener al menos dos estados")
    requested_state_cut = max(1, int(np.floor(n_neurons / 6.0 + 0.5)))
    state_cut = min(matrix.shape[0] - 1, requested_state_cut)
    if matrix.shape[0] > 1500:
        from scipy.sparse import csr_matrix
        from scipy.sparse.linalg import eigsh

        # M is symmetric. Its right singular vectors are its eigenvectors and
        # singular values are absolute eigenvalues, so Lanczos can compute only
        # the state_cut components that the source algorithm ever examines.
        eigenvalues, right_vectors = eigsh(
            csr_matrix((matrix + matrix.T) / 2.0),
            k=state_cut,
            which="LM",
            tol=1e-5,
            maxiter=max(1000, matrix.shape[0] * 2),
            v0=np.random.default_rng(0).normal(size=matrix.shape[0]),
        )
        singular_values = np.abs(eigenvalues)
        order = np.argsort(singular_values)[::-1]
        singular_values = singular_values[order]
        right_vectors = right_vectors[:, order]
        svd_method = "symmetric_eigensolver_svd_equivalent_top_state_cut"
    else:
        _, all_singular_values, all_vh = np.linalg.svd(matrix, full_matrices=False)
        singular_values = all_singular_values[:state_cut]
        vh = all_vh[:state_cut]
        right_vectors = vh.T
        svd_method = "dense_exact"
    factor_cut = float(initial_factor_cut)
    supports: list[np.ndarray] = []
    selected_factors = np.empty(0, dtype=int)
    loadings = np.zeros((0, matrix.shape[0]), dtype=float)

    for _ in range(1000):
        components = []
        counts = []
        for factor in range(state_cut):
            component = singular_values[factor] * np.outer(
                right_vectors[:, factor], right_vectors[:, factor]
            )
            selected = component > factor_cut
            components.append(selected)
            counts.append(int(selected.sum()))
        counts_array = np.asarray(counts, dtype=float)
        sizes = np.floor(np.sqrt(counts_array)).astype(int)
        total = int(sizes.sum())
        if total == 0:
            supports = []
            selected_factors = np.empty(0, dtype=int)
            loadings = np.zeros((0, matrix.shape[0]), dtype=float)
            break
        selected_factors = np.flatnonzero((sizes / total) >= minimum_prevalence)
        supports = [
            np.flatnonzero(np.any(components[factor], axis=0))
            for factor in selected_factors
        ]
        loadings = np.zeros((len(supports), matrix.shape[0]), dtype=float)
        for index, states in enumerate(supports):
            loadings[index, states] = np.abs(right_vectors[states, selected_factors[index]])
        if not supports:
            selected_factors = np.empty(0, dtype=int)
            break
        membership = np.zeros((len(supports), matrix.shape[0]), dtype=np.uint8)
        for index, states in enumerate(supports):
            membership[index, states] = 1
        if int(membership.sum(axis=0).max(initial=0)) <= 1:
            break
        factor_cut += factor_cut_step
    else:
        supports = []
        selected_factors = np.empty(0, dtype=int)
        loadings = np.zeros((0, matrix.shape[0]), dtype=float)

    # The paper's raster is ordered by first appearance of each state.
    order = np.argsort([states.min() if states.size else matrix.shape[0] for states in supports])
    supports = [supports[index] for index in order]
    loadings = loadings[order] if len(order) else loadings
    singular_ranks = selected_factors[order] + 1 if len(order) else np.empty(0, dtype=int)
    if return_singular_ranks:
        return singular_values, right_vectors, supports, factor_cut, state_cut, svd_method, singular_ranks
    return singular_values, right_vectors, supports, factor_cut, state_cut, svd_method


def _roc_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    labels = np.asarray(labels, dtype=bool)
    scores = np.asarray(scores, dtype=float)
    positives = int(labels.sum())
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        return float("nan")
    ranks = rankdata(scores, method="average")
    return float((ranks[labels].sum() - positives * (positives + 1) / 2) / (positives * negatives))


def identify_core_neurons(
    binary_population_vectors: np.ndarray,
    tfidf_population_vectors: np.ndarray,
    state_labels: np.ndarray,
    factor_supports: list[np.ndarray],
    neuron_ids: np.ndarray,
) -> list[dict[str, object]]:
    """Select representative cells by Stoixeion's in-sample ROC criterion.

    Labels and candidate cells come from the same recording. The separate
    within-day phase projection tests generalization; this AUC by itself is
    not independent validation of an ensemble.
    """
    results = []
    thresholds = np.arange(0.01, 0.101, 0.01)
    for factor, _states in enumerate(factor_supports, start=1):
        in_factor = state_labels == factor
        if not np.any(in_factor):
            results.append(
                {"factor": factor, "core_ids": [], "pool_ids": [], "auc": float("nan"), "threshold": float("nan")}
            )
            continue
        state_hist = tfidf_population_vectors[:, in_factor].sum(axis=1)
        maximum = float(state_hist.max(initial=0.0))
        if maximum <= 0:
            results.append(
                {"factor": factor, "core_ids": [], "pool_ids": [], "auc": float("nan"), "threshold": float("nan")}
            )
            continue
        state_hist = state_hist / maximum
        labels = in_factor
        best_auc = -np.inf
        best_threshold = float("nan")
        best_core = np.empty(0, dtype=int)
        for threshold in thresholds:
            core = np.flatnonzero(state_hist > threshold)
            if core.size == 0:
                continue
            query = np.zeros(binary_population_vectors.shape[0], dtype=float)
            query[core] = 1.0
            query_norm = np.linalg.norm(query)
            frame_vectors = binary_population_vectors.astype(float, copy=False)
            norms = np.linalg.norm(frame_vectors, axis=0)
            similarity = np.divide(
                query @ frame_vectors,
                query_norm * norms,
                out=np.zeros(frame_vectors.shape[1], dtype=float),
                where=norms > 0,
            )
            auc = _roc_auc(labels, similarity)
            if np.isfinite(auc) and auc > best_auc:
                best_auc = auc
                best_threshold = float(threshold)
                best_core = core
        pool = np.flatnonzero(state_hist > 0)
        results.append(
            {
                "factor": factor,
                "core_ids": [int(neuron_ids[index]) for index in best_core],
                "pool_ids": [int(neuron_ids[index]) for index in pool],
                "core_columns": best_core,
                "pool_columns": pool,
                "auc": float(best_auc) if np.isfinite(best_auc) else float("nan"),
                "threshold": best_threshold,
                "state_hist": state_hist,
            }
        )
    return results


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_condition(
    fluorescence: np.ndarray | None,
    deconvolved: np.ndarray | None,
    timestamps: np.ndarray | None,
    neuron_ids: np.ndarray,
    output_dir: Path,
    *,
    threshold_sd: float = 3.0,
    pks_shuffles: int = 100,
    pks_null_vector_cap: int = 2000,
    pks_null_percentile: float = 98.0,
    pks_cutoff_method: str = "continuous_percentile",
    scut: float | None = None,
    scut_shuffles: int = 100,
    hcut: float = 0.28,
    pks_override: int | None = None,
    pks_source: str = "automatic",
    seed: int = 1729,
    condition: str = "condition",
    events_override: np.ndarray | None = None,
    event_thresholds_override: np.ndarray | None = None,
    event_definition: str | None = None,
    save_matrices: bool = True,
    plot_phase_ranges: list[tuple[int, int]] | None = None,
    plot_phase_names: list[str] | None = None,
    plot_phase_durations_seconds: list[float] | None = None,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    if events_override is None:
        if fluorescence is None or deconvolved is None or timestamps is None:
            raise ValueError("C, S y t son obligatorios si no se entrega events_override")
        events, activity_thresholds, bin_info = binarize_from_calcium_derivative(
            fluorescence, deconvolved, timestamps, threshold_sd=threshold_sd
        )
        if event_definition is None:
            event_definition = "positive dC/dt above per-cell quiescent-frame mean + threshold SD"
    else:
        events = np.asarray(events_override)
        if events.ndim != 2 or events.shape[1] != len(neuron_ids):
            raise ValueError("events_override debe ser frames x neuronas alineadas con neuron_ids")
        events = (events > 0).astype(np.uint8)
        if event_thresholds_override is None:
            activity_thresholds = np.full(events.shape[1], np.nan, dtype=float)
        else:
            activity_thresholds = np.asarray(event_thresholds_override, dtype=float).reshape(-1)
            if len(activity_thresholds) != events.shape[1]:
                raise ValueError("event_thresholds_override no coincide con neuron_ids")
        finite_thresholds = activity_thresholds[np.isfinite(activity_thresholds)]
        median_threshold = float(np.median(finite_thresholds)) if finite_thresholds.size else float("nan")
        bin_info = BinarizationInfo(
            threshold_sd=float(threshold_sd),
            active_fraction=float(events.mean()),
            median_threshold=median_threshold,
            median_noise_mean=float("nan"),
            median_noise_sd=float("nan"),
            cells_using_fallback=0,
        )
        if event_definition is None:
            event_definition = "provided binary event raster"
    if pks_override is None:
        significant_frames, pks, pks_info = select_significant_vectors(
            events,
            n_shuffles=pks_shuffles,
            null_percentile=pks_null_percentile,
            seed=seed,
            max_null_vectors_per_shuffle=(
                None if pks_null_vector_cap <= 0 else pks_null_vector_cap
            ),
            cutoff_method=pks_cutoff_method,
        )
        pks_mode = "automatic_shuffle_on_control"
    else:
        pks = int(pks_override)
        significant_frames, pks_info = select_vectors_at_fixed_pks(events, pks, pks_source)
        pks_mode = "frozen_from_control"
    if significant_frames.size < 3:
        summary = {
            "condition": condition,
            "event_definition": event_definition,
            "status": "insufficient_significant_vectors",
            "threshold_sd": float(threshold_sd),
            "seed": int(seed),
            "n_neurons": int(events.shape[1]),
            "n_frames": int(events.shape[0]),
            "event_fraction": bin_info.active_fraction,
            "pks": int(pks),
            "pks_status": pks_info.get("status", "unknown"),
            "pks_shuffles": int(pks_shuffles),
            "pks_null_vector_cap": int(pks_null_vector_cap),
            "pks_null_percentile": float(pks_null_percentile),
            "pks_cutoff_method": pks_cutoff_method,
            "n_significant_vectors": int(significant_frames.size),
            "n_ensembles": 0,
        }
        (output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        (output_dir / "pks_diagnostics.json").write_text(json.dumps(pks_info, indent=2), encoding="utf-8")
        return {
            **summary,
            "cores": [],
            "neuron_ids": neuron_ids,
            "factor_rows": [],
            "significant_frames": significant_frames,
            "state_labels": np.zeros(len(significant_frames), dtype=int),
            "singular_values": np.empty(0, dtype=float),
        }

    binary_vectors = events[significant_frames].T.astype(np.uint8)
    tfidf_vectors, idf = stoixeion_tfidf(binary_vectors)
    similarity = cosine_similarity_map(tfidf_vectors)
    if scut is None:
        scut = estimate_scut_from_shuffles(
            tfidf_vectors,
            n_shuffles=scut_shuffles,
            percentile=88.0,
            seed=seed + 101,
        )
    initial_binary = (similarity > scut).astype(np.uint8)
    np.fill_diagonal(initial_binary, 1)
    del similarity
    initial_edges = int((initial_binary.sum() - len(initial_binary)) // 2)
    final_binary = two_pass_hamming_denoise(initial_binary, hcut=hcut)
    final_edges = int((final_binary.sum() - len(final_binary)) // 2)
    singular_values, right_vectors, supports, factor_cut, state_cut, svd_method, singular_ranks = svd_ensemble_support(
        final_binary, n_neurons=events.shape[1], return_singular_ranks=True
    )
    state_labels = np.zeros(len(significant_frames), dtype=int)
    for factor, states in enumerate(supports, start=1):
        state_labels[states] = factor
    tfidf_cores = identify_core_neurons(
        binary_vectors,
        tfidf_vectors,
        state_labels,
        supports,
        np.asarray(neuron_ids),
    )

    factor_rows = []
    core_rows = []
    for item in tfidf_cores:
        factor = int(item["factor"])
        support = supports[factor - 1]
        core_ids = item.get("core_ids", [])
        factor_rows.append(
            {
                "condition": condition,
                "factor": factor,
                "singular_rank": int(singular_ranks[factor - 1]),
                "singular_value": float(singular_values[singular_ranks[factor - 1] - 1]),
                "support_vectors": int(len(support)),
                "support_fraction_of_significant_vectors": float(len(support) / len(significant_frames)),
                "support_fraction_of_recording": float(len(support) / len(events)),
                "passes_5pct_full_recording_sensitivity": bool(len(support) / len(events) >= 0.05),
                "pool_size": int(len(item.get("pool_ids", []))),
                "core_size": int(len(core_ids)),
                "core_auc": float(item.get("auc", float("nan"))),
                "core_threshold": float(item.get("threshold", float("nan"))),
            }
        )
        for rank, neuron_id in enumerate(core_ids, start=1):
            core_rows.append(
                {
                    "condition": condition,
                    "factor": factor,
                    "rank": rank,
                    "neuron_id_within_merged_day": int(neuron_id),
                }
            )

    summary = {
        "condition": condition,
        "event_definition": event_definition,
        "status": "ok" if supports else "no_ensemble_passed_svd_cutoffs",
        "threshold_sd": float(threshold_sd),
        "seed": int(seed),
        "n_neurons": int(events.shape[1]),
        "n_frames": int(events.shape[0]),
        "event_fraction": bin_info.active_fraction,
        "median_event_threshold": bin_info.median_threshold,
        "median_noise_mean": bin_info.median_noise_mean,
        "median_noise_sd": bin_info.median_noise_sd,
        "binarization_fallback_cells": bin_info.cells_using_fallback,
        "pks": int(pks),
        "pks_selection_mode": pks_mode,
        "pks_source": pks_source,
        "pks_null_percentile": float(pks_null_percentile),
        "pks_cutoff_method": pks_cutoff_method,
        "pks_shuffles": int(pks_shuffles),
        "pks_null_vector_cap": int(pks_null_vector_cap),
        "pks_status": pks_info.get("status", "unknown"),
        "pks_real_mean_similarity": pks_info.get("real_mean_similarity", float("nan")),
        "pks_shuffle_cutoff": pks_info.get("shuffle_similarity_cutoff", float("nan")),
        "n_significant_vectors": int(len(significant_frames)),
        "tfidf": True,
        "scut": float(scut),
        "scut_method": "automatic Stoixeion.calc_scut; 88th percentile of time-shuffled cosine similarities",
        "scut_shuffles": int(scut_shuffles),
        "scut_pair_sampling_used": bool(binary_vectors.shape[1] > 1500),
        "scut_pair_samples_per_shuffle": 100_000 if binary_vectors.shape[1] > 1500 else 0,
        "hcut": float(hcut),
        "n_high_similarity_edges_before_denoising": initial_edges,
        "n_final_similarity_edges": final_edges,
        "state_cut_upper_bound": int(state_cut),
        "svd_method": svd_method,
        "factor_cut_final": float(factor_cut),
        "n_ensembles": int(len(supports)),
        "n_ensembles_5pct_full_recording_sensitivity": int(sum(
            len(support) / len(events) >= 0.05 for support in supports
        )),
        "singular_values_top10": [float(x) for x in singular_values[:10]],
        "factor_singular_ranks": [int(value) for value in singular_ranks],
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (output_dir / "pks_diagnostics.json").write_text(json.dumps(pks_info, indent=2), encoding="utf-8")
    if save_matrices:
        np.savez_compressed(
            output_dir / "carrillo_reid_matrices.npz",
            binary_activity=events,
            significant_frames=significant_frames,
            binary_population_vectors=binary_vectors,
            tfidf_population_vectors=tfidf_vectors.astype(np.float32),
            idf=idf.astype(np.float32),
            similarity_binary=initial_binary,
            similarity_final=final_binary,
            singular_values=singular_values,
            state_labels=state_labels,
            neuron_ids=np.asarray(neuron_ids),
            event_thresholds=activity_thresholds,
        )
    _write_csv(
        output_dir / "ensemble_summary.csv",
        list(factor_rows[0]) if factor_rows else [
            "condition", "factor", "singular_rank", "singular_value", "support_vectors",
            "support_fraction_of_significant_vectors", "support_fraction_of_recording",
            "passes_5pct_full_recording_sensitivity",
            "pool_size", "core_size", "core_auc", "core_threshold",
        ],
        factor_rows,
    )
    _write_csv(
        output_dir / "core_neurons.csv",
        ["condition", "factor", "rank", "neuron_id_within_merged_day"],
        core_rows,
    )
    _plot_condition_clean(
        events,
        initial_binary,
        singular_values,
        state_labels,
        condition,
        output_dir / "carrillo_reid_figure.png",
        significant_frames=significant_frames,
        phase_ranges=plot_phase_ranges,
        phase_names=plot_phase_names,
        phase_durations_seconds=plot_phase_durations_seconds,
    )
    if plot_phase_ranges and plot_phase_names:
        _plot_condition_clean(
            events,
            initial_binary,
            singular_values,
            state_labels,
            condition,
            output_dir / "carrillo_reid_figure_vector_order.png",
            significant_frames=significant_frames,
            phase_ranges=plot_phase_ranges,
            phase_names=plot_phase_names,
            phase_durations_seconds=plot_phase_durations_seconds,
            d_axis_mode="vector_order",
        )
    return {
        **summary,
        "cores": [set(item.get("core_ids", [])) for item in tfidf_cores],
        "neuron_ids": np.asarray(neuron_ids),
        "factor_rows": factor_rows,
        "significant_frames": significant_frames,
        "state_labels": state_labels,
        "singular_values": singular_values,
    }


def _plot_condition_clean(
    events: np.ndarray,
    initial_binary: np.ndarray,
    singular_values: np.ndarray,
    state_labels: np.ndarray,
    condition: str,
    output_path: Path,
    *,
    significant_frames: np.ndarray | None = None,
    phase_ranges: list[tuple[int, int]] | None = None,
    phase_names: list[str] | None = None,
    phase_durations_seconds: list[float] | None = None,
    d_axis_mode: str = "frame",
) -> None:
    def display_map(binary_map: np.ndarray, max_side: int = 1200) -> np.ndarray:
        n = binary_map.shape[0]
        if n <= max_side:
            return binary_map
        block = int(np.ceil(n / max_side))
        pad = (-n) % block
        padded = np.pad(binary_map, ((0, pad), (0, pad)), mode="constant")
        return padded.reshape(
            padded.shape[0] // block,
            block,
            padded.shape[1] // block,
            block,
        ).sum(axis=(1, 3), dtype=np.float32) / (block * block)

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    axes[0, 0].imshow(events.T, aspect="auto", interpolation="nearest", cmap="Greys", vmin=0, vmax=1)
    axes[0, 0].set_title("A: actividad binaria")
    axes[0, 0].set_xlabel("Frame")
    axes[0, 0].set_ylabel("Neurona")

    if phase_ranges and phase_names:
        phase_colors = plt.get_cmap("Set2")
        for phase_index, ((start, stop), phase_name) in enumerate(
            zip(phase_ranges, phase_names)
        ):
            color = phase_colors(phase_index % phase_colors.N)
            left = max(float(start) - 0.5, -0.5)
            right = min(float(stop) - 0.5, events.shape[0] - 0.5)
            if right <= left:
                continue
            axes[0, 0].axvspan(left, right, color=color, alpha=0.10, zorder=0.5)
            if phase_index < len(phase_ranges) - 1:
                axes[0, 0].axvline(right, color="black", linewidth=0.8, alpha=0.65, zorder=3)

    axes[0, 1].imshow(display_map(initial_binary), aspect="auto", interpolation="nearest", cmap="Greys", vmin=0, vmax=1)
    axes[0, 1].set_title("B: M binaria")
    axes[0, 1].set_xlabel("Vector poblacional significativo")
    axes[0, 1].set_ylabel("Vector poblacional significativo")

    ranks = np.arange(1, len(singular_values) + 1)
    axes[1, 0].plot(ranks, singular_values, color="black", linewidth=1)
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_title("C: valores singulares")
    axes[1, 0].set_xlabel("Rango singular")
    axes[1, 0].set_ylabel("Valor singular")

    have_frame_map = (
        significant_frames is not None
        and len(significant_frames) == len(state_labels)
    )
    frame_positions = (
        np.asarray(significant_frames, dtype=float)
        if have_frame_map
        else np.arange(len(state_labels), dtype=float)
    )
    if d_axis_mode not in {"frame", "vector_order"}:
        raise ValueError("d_axis_mode must be 'frame' or 'vector_order'")
    if d_axis_mode == "frame":
        d_positions = frame_positions
        d_phase_ranges = phase_ranges or []
        d_xlabel = "Frame concatenado" if have_frame_map else "Vector significativo (orden temporal)"
    else:
        d_positions = np.arange(len(state_labels), dtype=float)
        d_xlabel = "Vector poblacional significativo (orden temporal)"
        d_phase_ranges = []
        if phase_ranges and have_frame_map:
            cursor = 0
            for start, stop in phase_ranges:
                n_vectors = int(np.sum((frame_positions >= start) & (frame_positions < stop)))
                d_phase_ranges.append((cursor, cursor + n_vectors))
                cursor += n_vectors

    if np.any(state_labels):
        n_factors = int(state_labels.max())
        if d_axis_mode == "vector_order":
            raster = np.zeros((n_factors, len(state_labels)), dtype=np.uint8)
            for factor in range(1, n_factors + 1):
                raster[factor - 1, state_labels == factor] = 1
            axes[1, 1].imshow(
                raster,
                aspect="auto",
                interpolation="nearest",
                cmap=ListedColormap([(0, 0, 0, 0), (0, 0, 0, 1)]),
                vmin=0,
                vmax=1,
                zorder=1,
            )
        else:
            for factor in range(1, n_factors + 1):
                event_frames = d_positions[state_labels == factor]
                row = factor - 1
                axes[1, 1].vlines(
                    event_frames,
                    row - 0.32,
                    row + 0.32,
                    color="black",
                    linewidth=0.45,
                    zorder=2,
                )
        axes[1, 1].set_yticks(np.arange(n_factors), labels=[f"E{i}" for i in range(1, n_factors + 1)])
        axes[1, 1].set_ylim(n_factors - 0.5, -0.5)
        axes[1, 1].set_xlabel(d_xlabel)
        if d_positions.size:
            if d_axis_mode == "frame" and have_frame_map:
                axes[1, 1].set_xlim(-0.5, events.shape[0] - 0.5)
            else:
                axes[1, 1].set_xlim(-0.5, len(state_labels) - 0.5)
    elif d_axis_mode == "frame" and have_frame_map:
        axes[1, 1].set_xlim(-0.5, events.shape[0] - 0.5)

    if d_phase_ranges and phase_names:
        phase_centers = []
        phase_tick_labels = []
        phase_colors = plt.get_cmap("Set2")
        for phase_index, ((start, stop), phase_name) in enumerate(
            zip(d_phase_ranges, phase_names)
        ):
            if stop <= start:
                continue
            color = phase_colors(phase_index % phase_colors.N)
            axes[1, 1].axvspan(
                start - 0.5,
                stop - 0.5,
                color=color,
                alpha=0.10,
                zorder=0.5,
            )
            if phase_index < len(d_phase_ranges) - 1:
                axes[1, 1].axvline(
                    stop - 0.5,
                    color="black",
                    linewidth=0.8,
                    alpha=0.65,
                    zorder=3,
                )
            phase_centers.append((start + stop - 1) / 2.0)
            label = phase_name
            if phase_durations_seconds and phase_index < len(phase_durations_seconds):
                duration_min = phase_durations_seconds[phase_index] / 60.0
                label = f"{phase_name}\n{duration_min:.1f} min"
            phase_tick_labels.append(label)
        if phase_centers:
            phase_axis = axes[1, 1].twiny()
            phase_axis.set_xlim(axes[1, 1].get_xlim())
            phase_axis.set_xticks(phase_centers)
            phase_axis.set_xticklabels(phase_tick_labels, fontsize=8)
            phase_axis.tick_params(axis="x", length=0, pad=2)
            phase_axis.set_xlabel("Subsesión", labelpad=2)

    axes[1, 1].set_title("D: vectores asignados a ensembles")
    axes[1, 1].set_ylabel("Ensemble")

    fig.suptitle(f"Carrillo-Reid / Stoixeion | {condition}", fontsize=14)
    fig.savefig(output_path, dpi=170)
    plt.close(fig)


def _plot_condition(
    events: np.ndarray,
    significant_frames: np.ndarray,
    initial_binary: np.ndarray,
    final_binary: np.ndarray,
    singular_values: np.ndarray,
    state_labels: np.ndarray,
    condition: str,
    scut: float,
    hcut: float,
    output_path: Path,
) -> None:
    def display_map(binary_map: np.ndarray, max_side: int = 1200) -> tuple[np.ndarray, int]:
        """Block-average a binary map for plotting only; never for analysis."""
        n = binary_map.shape[0]
        if n <= max_side:
            return binary_map, 1
        block = int(np.ceil(n / max_side))
        pad = (-n) % block
        padded = np.pad(binary_map, ((0, pad), (0, pad)), mode="constant")
        reduced = padded.reshape(padded.shape[0] // block, block,
                                 padded.shape[1] // block, block).sum(
                                     axis=(1, 3), dtype=np.float32
                                 ) / (block * block)
        return reduced, block

    initial_display, initial_block = display_map(initial_binary)
    final_display, final_block = display_map(final_binary)
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    axes[0, 0].imshow(events.T, aspect="auto", interpolation="nearest", cmap="Greys", vmin=0, vmax=1)
    axes[0, 0].set_title("A: actividad binaria (dC/dt > basal + 3 SD)")
    axes[0, 0].set_xlabel("Frame")
    axes[0, 0].set_ylabel("Neurona")

    axes[0, 1].imshow(initial_display, aspect="auto", interpolation="nearest", cmap="Greys", vmin=0, vmax=1)
    axes[0, 1].set_title(f"B: M binaria (scut={scut:.2f}; visual {initial_block}×)")
    axes[0, 1].set_xlabel("Vector poblacional significativo")
    axes[0, 1].set_ylabel("Vector poblacional significativo")

    axes[0, 2].imshow(final_display, aspect="auto", interpolation="nearest", cmap="Greys", vmin=0, vmax=1)
    axes[0, 2].set_title(f"C: M final (hcut={hcut:.2f}; visual {final_block}×)")
    axes[0, 2].set_xlabel("Vector poblacional")
    axes[0, 2].set_ylabel("Vector poblacional")

    ranks = np.arange(1, len(singular_values) + 1)
    axes[1, 0].plot(ranks, singular_values, color="black", linewidth=1)
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_title("D: valores singulares usados por Stoixeion")
    axes[1, 0].set_xlabel("Rango singular")
    axes[1, 0].set_ylabel("Valor singular (escala log)")

    if np.any(state_labels):
        raster = np.zeros((int(state_labels.max()), len(state_labels)), dtype=np.uint8)
        for factor in range(1, int(state_labels.max()) + 1):
            raster[factor - 1, state_labels == factor] = 1
            axes[1, 1].imshow(raster, aspect="auto", interpolation="nearest", cmap="Greys", vmin=0, vmax=1)
    axes[1, 1].set_title("E: vectores asignados a cada ensemble")
    axes[1, 1].set_xlabel("Vector poblacional significativo (orden temporal)")
    axes[1, 1].set_ylabel("Ensemble")
    axes[1, 2].axis("off")
    axes[1, 2].text(
        0.02,
        0.95,
        "F: parámetros\n"
        f"frames: {len(events):,}\n"
        f"neuronas: {events.shape[1]}\n"
        f"vectores significativos: {len(significant_frames):,}\n"
        f"eventos binarios: {events.mean():.1%}\n"
        f"ensembles: {int(state_labels.max(initial=0))}\n"
        "Mapas en orden temporal\n"
        "El agrupamiento por bloques afecta solo la figura",
        transform=axes[1, 2].transAxes,
        va="top",
        ha="left",
        fontsize=12,
        color="black",
    )
    fig.suptitle(f"Carrillo-Reid / Stoixeion | {condition} | {len(significant_frames)} vectores", fontsize=14)
    fig.savefig(output_path, dpi=170)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Carrillo-Reid Stoixeion pipeline on one merged subsession.")
    parser.add_argument("--animal", required=True)
    parser.add_argument("--session", required=True)
    parser.add_argument("--subsession", type=int, required=True)
    parser.add_argument("--condition", default=None)
    parser.add_argument("--threshold-sd", type=float, default=3.0)
    parser.add_argument("--pks-shuffles", type=int, default=100)
    parser.add_argument("--pks-cutoff-method", choices=("continuous_percentile", "stoixeion_histc"), default="continuous_percentile")
    parser.add_argument("--pks-null-vector-cap", type=int, default=2000,
                        help="0 uses all null vectors, as in the MATLAB source")
    parser.add_argument("--scut", default="auto", help="auto (Stoixeion) or a numeric cutoff")
    parser.add_argument("--scut-shuffles", type=int, default=100)
    parser.add_argument("--hcut", type=float, default=0.28)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--output-dir", default="coactivation_analysis/results/current/carrillo_reid")
    args = parser.parse_args()

    from core.data_loader import get_common_neurons, load_session, mapping_columns

    data = load_session(args.animal, args.session)
    if data["type"] != "merged":
        raise ValueError("Carrillo-Reid batch currently requires merged four-condition files")
    sub_index = args.subsession - 1
    sub = data["subsessions"][sub_index]
    if sub.get("C") is None:
        raise ValueError("act.C is required for the 3 SD fluorescence-derivative threshold")
    neuron_ids = get_common_neurons(data["mapping"])
    columns = mapping_columns(data["mapping"], neuron_ids, sub_index, sub["S"].shape[1])
    output_dir = Path(args.output_dir) / f"{args.animal}_{Path(args.session).stem}" / f"subsession_{args.subsession}"
    condition = args.condition or f"subsession_{args.subsession}"
    scut = None if args.scut.lower() == "auto" else float(args.scut)
    result = run_condition(
        sub["C"][:, columns],
        sub["S"][:, columns],
        sub["t"],
        neuron_ids,
        output_dir,
        threshold_sd=args.threshold_sd,
        pks_shuffles=args.pks_shuffles,
        pks_cutoff_method=args.pks_cutoff_method,
        pks_null_vector_cap=args.pks_null_vector_cap,
        scut=scut,
        scut_shuffles=args.scut_shuffles,
        hcut=args.hcut,
        seed=args.seed,
        condition=condition,
    )
    print(json.dumps({key: value for key, value in result.items() if key not in {"cores", "neuron_ids", "factor_rows"}}, indent=2))
    print(f"Saved: {output_dir}")


if __name__ == "__main__":
    main()
