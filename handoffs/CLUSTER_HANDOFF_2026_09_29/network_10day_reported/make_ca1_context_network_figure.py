"""Poster-ready descriptive joint figure: spatial maps and CA1 networks.

The two axes answer different questions: a cross-day HabL reference for the
weakest map pair, and within-day SAMPLE-to-TEST network identity against nulls.
This is not a drug effect or a behavior-neural correlation.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "resultados_correlacion"
MAP_SOURCE = (
    ROOT / "results" / "place_cell_stability_four_phases" /
    "sensitivity_cross200" / "strict_fourphase_six_pair_null_margin_by_day.csv"
)
NETWORK_SOURCE = OUT / "All_object_networks" / "all_object_networks_onsets.csv"
NETWORK_FRAME_SOURCE = OUT / "All_object_networks" / "all_object_networks_frames.csv"
POSITION_SOURCE = OUT / "All_object_networks" / "position_controlled_networks_onsets.csv"
POSITION_FRAME_SOURCE = OUT / "All_object_networks" / "position_controlled_networks_frames.csv"


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def build_rows() -> list[dict[str, object]]:
    maps = _read_csv(MAP_SOURCE)
    networks = _read_csv(NETWORK_SOURCE)
    network_frames = {r["day_id"]: r for r in _read_csv(NETWORK_FRAME_SOURCE)}
    controlled = {r["day_id"]: r for r in _read_csv(POSITION_SOURCE)}
    controlled_frames = {r["day_id"]: r for r in _read_csv(POSITION_FRAME_SOURCE)}
    if len(networks) != 10:
        raise ValueError("expected ten object days with valid four-phase mapping")
    baseline = {
        (r["animal"], r["representation"]): float(r["median_min_r_minus_null_q95"])
        for r in maps if r["day_type"] == "baseline"
    }
    by_day = {(r["day_id"], r["representation"]): r for r in maps
              if r["day_type"] == "objects"}
    rows = []
    for net in networks:
        day_id = net["day_id"]
        animal = net["animal"]
        metrics = {}
        for representation, label in (
            ("binary S>3 SD maps", "event_map"),
            ("continuous max(S,0) maps", "continuous_map"),
        ):
            current = float(by_day[(day_id, representation)]["median_min_r_minus_null_q95"])
            metrics[f"delta_{label}_weakest_margin_vs_HabL"] = (
                current - baseline[(animal, representation)]
            )
        observed = float(net["network_similarity_observed"])
        null_max = max(float(net["identity_null_q95"]), float(net["temporal_null_q95"]))
        frame = network_frames[day_id]
        residual = controlled[day_id]
        residual_frame = controlled_frames[day_id]
        residual_null_max = max(float(residual["identity_null_q95"]),
                                float(residual["temporal_null_q95"]))
        rows.append({
            "animal": animal, "task": net["task"], "treatment": net["treatment"],
            "day_id": day_id,
            **metrics,
            "network_onset_similarity_observed": observed,
            "network_onset_max_null_q95": null_max,
            "network_onset_observed_over_max_null_q95": observed / null_max,
            "network_onset_identity_holm_p": float(net["identity_holm_p"]),
            "network_onset_temporal_holm_p": float(net["temporal_holm_p"]),
            "network_frame_similarity_observed": float(frame["network_similarity_observed"]),
            "network_frame_max_null_q95": max(
                float(frame["identity_null_q95"]), float(frame["temporal_null_q95"])
            ),
            "position_controlled_network_onset_similarity": float(
                residual["network_similarity_observed"]
            ),
            "position_controlled_network_onset_max_null_q95": residual_null_max,
            "position_controlled_network_onset_ratio": (
                float(residual["network_similarity_observed"]) / residual_null_max
            ),
            "position_controlled_identity_holm_p": float(residual["identity_holm_p"]),
            "position_controlled_temporal_holm_p": float(residual["temporal_holm_p"]),
            "position_controlled_network_frame_similarity": float(
                residual_frame["network_similarity_observed"]
            ),
            "position_controlled_network_frame_max_null_q95": max(
                float(residual_frame["identity_null_q95"]),
                float(residual_frame["temporal_null_q95"]),
            ),
        })
    if len({r["day_id"] for r in rows}) != 10:
        raise AssertionError("duplicate or missing object day")
    if not all(
        r["delta_event_map_weakest_margin_vs_HabL"] < 0
        and r["delta_continuous_map_weakest_margin_vs_HabL"] < 0
        and r["network_onset_similarity_observed"] > r["network_onset_max_null_q95"]
        and r["network_frame_similarity_observed"] > r["network_frame_max_null_q95"]
        and r["position_controlled_network_onset_similarity"]
            > r["position_controlled_network_onset_max_null_q95"]
        and r["position_controlled_network_frame_similarity"]
            > r["position_controlled_network_frame_max_null_q95"]
        for r in rows
    ):
        raise AssertionError("not all ten days reproduce both descriptive findings")
    return rows


def plot(rows: list[dict[str, object]], path: Path) -> None:
    rows = sorted(rows, key=lambda r: (
        r["animal"], 0 if r["task"] == "SD" else 1,
        0 if r["treatment"] == "VEH" else 1,
    ))
    y = np.arange(len(rows))[::-1]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 7), sharey=True,
                             gridspec_kw={"width_ratios": [1, 1]})
    colors = {"CNO": "#326aa8", "VEH": "#db7d28"}
    symbols = {"SD": "o", "XsS": "s"}
    for pos, row in zip(y, rows):
        color = colors[row["treatment"]]
        marker = symbols[row["task"]]
        axes[0].scatter(row["delta_event_map_weakest_margin_vs_HabL"], pos,
                        color=color, marker=marker, s=90, zorder=3)
        axes[1].scatter(row["position_controlled_network_onset_ratio"], pos,
                        color=color, marker=marker, s=90, zorder=3)
    labels = [f"{r['animal']}  {r['task']}  {r['treatment']}" for r in rows]
    axes[0].set_yticks(y, labels)
    axes[0].axvline(0, color="#555555", linestyle="--", linewidth=1)
    axes[1].axvline(1, color="#555555", linestyle="--", linewidth=1)
    axes[0].set_xlim(-.15, .015)
    axes[1].set_xlim(0, max(r["network_onset_observed_over_max_null_q95"] for r in rows) * 1.12)
    axes[0].set_xlabel("Cambio del margen del mapa más débil vs HabL")
    axes[1].set_xlabel("Red tras controlar posición/velocidad / corte nulo 95%")
    axes[0].set_title("Mapas menos repetibles: 10/10 jornadas")
    axes[1].set_title("Coactivación residual: 10/10 jornadas")
    for ax in axes:
        ax.grid(axis="x", alpha=.2)
        ax.set_ylim(-.7, len(rows) - .3)
    fig.suptitle("CA1 durante tareas con objetos: cambian los mapas, persiste la coactividad",
                 y=.96)
    fig.subplots_adjust(left=.24, right=.98, top=.86, bottom=.18, wspace=.08)
    fig.text(.61, .05,
             "Tres animales, diez jornadas; el baseline HabL es otro día. "
             "Red: 199 nulos por tipo, p Holm=0,05 entre diez jornadas; frames replican. "
             "No es efecto CNO ni correlación conductual individual.",
             ha="center", fontsize=9, wrap=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    rows = build_rows()
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "CA1_context_network_joint_evidence.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    plot(rows, OUT / "18_CA1_maps_vs_networks.png")
    print("Verified map and network pattern in", len(rows), "object days.")
    print("Minimum onset network/null ratio:", round(min(
        r["position_controlled_network_onset_ratio"] for r in rows
    ), 2))
    print("Maximum Holm-adjusted p:", max(
        max(r["network_onset_identity_holm_p"],r["network_onset_temporal_holm_p"])
        for r in rows
    ))


if __name__ == "__main__":
    main()
