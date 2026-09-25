import argparse
import os
import re
from pathlib import Path

import matplotlib
import numpy as np
import scipy.io

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from core.data_loader import (
    get_common_neurons,
    get_data_dir,
    list_sessions,
    load_session,
    mapping_columns,
)
from core.place_fields import gaussian_rate_maps_matlab


ARENA_RADIUS_CM = 25.0
BIN_SIZE_CM = 1.0
SIGMA_CM = 4.0


def _phase_names(day_name, n_subsessions):
    if n_subsessions == 1:
        return [day_name or "Registro"]
    if day_name == "HabL":
        return ["OF1", "OF2", "OF3", "OF4"]
    if "_SD_" in day_name or "_XsS_" in day_name:
        return ["OF1", "SAMPLE", "TEST", "OF2"]
    return [f"Sub {index + 1}" for index in range(n_subsessions)]


def _prepare_rate_maps(subsessions):
    """Use the same continuous-S Gaussian rate-map recipe as the HabL gallery."""
    caches = []
    for sub in subsessions:
        x = np.asarray(sub["x"], dtype=float)
        y = np.asarray(sub["y"], dtype=float)
        t = np.asarray(sub["t"], dtype=float)
        raw_s = np.asarray(sub["S"], dtype=float)
        valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(t)
        valid &= np.isfinite(raw_s).all(axis=1)
        x, y, t = x[valid], y[valid], t[valid]
        positive_s = np.maximum(raw_s[valid], 0.0).astype(np.float32)
        rates, edges, meta = gaussian_rate_maps_matlab(
            x, y, t, positive_s,
            bin_size=BIN_SIZE_CM,
            sigma=SIGMA_CM,
            arena_radius=ARENA_RADIUS_CM,
        )
        caches.append({
            "x": x, "y": y, "t": t, "S": positive_s,
            "rates": rates, "edges": edges,
            "sample_rate_hz": meta["sample_rate_hz"],
        })
    return caches


def _style_axis(ax):
    ax.set_aspect("equal")
    ax.set_xlim(-ARENA_RADIUS_CM, ARENA_RADIUS_CM)
    ax.set_ylim(-ARENA_RADIUS_CM, ARENA_RADIUS_CM)
    ax.set_xticks((-25, 0, 25))
    ax.set_yticks((-25, 0, 25))
    ax.tick_params(labelsize=6)


def _render_cell_figure(animal, cell_label, cell_filename, columns,
                        caches, phase_names, output_dir):
    """Draw trajectory/S amplitude and continuous-S rate map; save PNG only."""
    n_phases = len(caches)
    fig, axes = plt.subplots(2, n_phases, figsize=(3.2 * n_phases, 6.2), squeeze=False)
    present = []
    for phase_index, cache in enumerate(caches):
        ax_trajectory = axes[0, phase_index]
        ax_map = axes[1, phase_index]
        _style_axis(ax_trajectory)
        _style_axis(ax_map)
        local_column, local_id = columns[phase_index]
        signal = cache["S"][:, local_column]
        positive = signal > 0
        n_positive = int(positive.sum())
        rate_map = cache["rates"][local_column]
        finite_rate = rate_map[np.isfinite(rate_map)]
        peak_rate = float(np.max(finite_rate)) if finite_rate.size else 0.0
        present.append(phase_names[phase_index])

        ax_trajectory.plot(
            cache["x"], cache["y"], color="0.82", linewidth=.3,
            rasterized=True, zorder=1,
        )
        if n_positive:
            scatter = ax_trajectory.scatter(
                cache["x"][positive], cache["y"][positive],
                c=signal[positive], cmap="turbo", vmin=0,
                vmax=float(np.max(signal[positive])), s=4.5, alpha=.85,
                linewidths=0, rasterized=True, zorder=2,
            )
            cb = fig.colorbar(scatter, ax=ax_trajectory, fraction=.046, pad=.03)
            cb.ax.tick_params(labelsize=6)
            cb.set_label("S", fontsize=6)
        ax_trajectory.set_title(
            f"{phase_names[phase_index]} (local {local_id})\n{n_positive} S > 0 samples",
            fontsize=8,
        )

        edges = cache["edges"]
        cmap = plt.cm.turbo.copy()
        cmap.set_bad(color="white")
        image = ax_map.imshow(
            rate_map.T, origin="lower",
            extent=[edges[0], edges[-1], edges[0], edges[-1]],
            cmap=cmap, vmin=0, vmax=peak_rate if peak_rate > 0 else 1,
            interpolation="nearest", rasterized=True,
        )
        cb = fig.colorbar(image, ax=ax_map, fraction=.046, pad=.03)
        cb.ax.tick_params(labelsize=6)
        cb.set_label("act/s", fontsize=6)
        ax_map.set_title(f"Activity rate map\nmax: {peak_rate:.2f} act/s", fontsize=8)

    fig.suptitle(
        f"{animal} | {cell_label} | Tracked in {n_phases}/{n_phases} phases "
        f"({','.join(present)})",
        fontsize=11, fontweight="semibold", y=.985,
    )
    fig.tight_layout(rect=(0, 0, 1, .94), w_pad=.8, h_pad=.8)
    fig.savefig(os.path.join(output_dir, cell_filename), dpi=130, facecolor="white")
    plt.close(fig)


def _process_local_only(animal, subsessions, session_out, phase_names):
    """Draw local cells per phase if cross-phase CellReg mapping is ambiguous."""
    total = 0
    for phase_index, (sub, cache) in enumerate(zip(subsessions, _prepare_rate_maps(subsessions))):
        local_out = os.path.join(session_out, f"phase_{phase_index + 1}_local_only")
        os.makedirs(local_out, exist_ok=True)
        columns = [(local_column, local_column) for local_column in range(sub["S"].shape[1])]
        for local_column in range(sub["S"].shape[1]):
            _render_cell_figure(
                animal, f"local column {local_column}",
                f"local_{local_column:04d}.png", [columns[local_column]],
                [cache], [phase_names[phase_index]], local_out,
            )
            total += 1
    print(
        f"  Guardados {total} mapas locales en {session_out}; "
        "no se enlazaron IDs entre fases.",
        flush=True,
    )


def process_session(animal, session_file, data_dir=None, out_dir="results",
                    category_parts=None):
    print(f"Procesando {animal} - {session_file}...", flush=True)
    data = load_session(animal, session_file, data_dir=data_dir)
    session_out = os.path.join(
        out_dir, animal, *(category_parts or ("Otros",)),
        session_file.replace(".mat", ""),
    )
    os.makedirs(session_out, exist_ok=True)
    subsessions = [data] if data["type"] == "simple" else data["subsessions"]
    day_name = str(getattr(data.get("sess"), "day_name", ""))
    phases = _phase_names(day_name, len(subsessions))

    if data["type"] == "simple":
        common_neurons = np.arange(data["n_neurons"], dtype=int)
        columns_per_sub = [common_neurons.copy()]
        mapping = None
    else:
        mapping = np.asarray(data["mapping"], dtype=float)
        common_neurons = get_common_neurons(mapping)
        try:
            columns_per_sub = [
                mapping_columns(mapping, common_neurons, index, sub["S"].shape[1])
                for index, sub in enumerate(subsessions)
            ]
        except ValueError as error:
            print(
                f"  CellReg no alinea todas las fases ({error}). "
                "Se guardarán mapas por fase, sin enlazar IDs.",
                flush=True,
            )
            _process_local_only(animal, subsessions, session_out, phases)
            return

    if len(common_neurons) == 0:
        if data["type"] == "merged":
            print("  Sin IDs comunes a todas las fases; mapas locales por fase.", flush=True)
            _process_local_only(animal, subsessions, session_out, phases)
            return
        print(f"  Sin neuronas: {animal} - {session_file}", flush=True)
        return

    caches = _prepare_rate_maps(subsessions)
    common_index = {int(global_id): index for index, global_id in enumerate(common_neurons)}
    for global_id in common_neurons:
        columns = []
        for phase_index in range(len(subsessions)):
            local_column = int(columns_per_sub[phase_index][common_index[int(global_id)]])
            local_id = (int(global_id) if mapping is None
                        else int(mapping[global_id, phase_index]))
            columns.append((local_column, local_id))
        _render_cell_figure(
            animal, f"Global cell {int(global_id):04d}",
            f"neuron_{int(global_id):04d}.png", columns, caches, phases,
            session_out,
        )

    print(f"Listos {len(common_neurons)} mapas en {session_out}", flush=True)


def _session_day_name(animal, session_file, data_dir):
    path = Path(data_dir) / animal / session_file
    metadata = scipy.io.loadmat(
        str(path), variable_names=["sess"], squeeze_me=True,
        struct_as_record=False,
    )
    sess = metadata.get("sess")
    return str(getattr(sess, "day_name", "")) if sess is not None else ""


def _session_category(day_name):
    """Simple output folders derived from sess.day_name."""
    if day_name == "HabL":
        return ("HabL",)
    if re.fullmatch(r"HabC\d*", day_name):
        return (day_name,)
    match = re.fullmatch(r"T\d+_(SD|XsS)_(CNO|VEH)_?", day_name)
    if match:
        task, treatment = match.groups()
        return ("Objetos", task, treatment)
    return ("Otros",)


def process_all(data_dir=None, out_dir="results/Mapas_por_tipo", animal=None, session_file=None,
                exclude_day_names=None):
    data_dir = data_dir or get_data_dir()
    excluded = set(exclude_day_names or ())
    sessions = [
        session for session in list_sessions(data_dir=data_dir)
        if (animal is None or session["animal"] == animal)
        and (session_file is None or session["session_file"] == session_file)
    ]
    categorized = []
    for session in sessions:
        try:
            day_name = _session_day_name(
                session["animal"], session["session_file"], data_dir
            )
        except (OSError, ValueError, KeyError):
            day_name = ""
        if day_name in excluded:
            print(
                f"[EXCLUIDA] {session['animal']} / {session['session_file']} "
                f"(day_name={day_name})",
                flush=True,
            )
            continue
        categorized.append((session, _session_category(day_name)))
    if not categorized:
        raise SystemExit("No se encontraron sesiones que coincidan con los filtros.")

    for session, category_parts in categorized:
        try:
            process_session(
                session["animal"],
                session["session_file"],
                data_dir=data_dir,
                out_dir=out_dir,
                category_parts=category_parts,
            )
        except Exception as error:
            print(
                f"[OMITIDA] {session['animal']} / {session['session_file']}: {error}",
                flush=True,
            )


def main():
    parser = argparse.ArgumentParser(
        description="Genera mapas de actividad para todas las sesiones disponibles."
    )
    parser.add_argument("--animal", help="Procesar solamente un animal")
    parser.add_argument("--session", dest="session_file", help="Procesar solamente un archivo .mat")
    parser.add_argument("--data-dir", default=None, help="Carpeta con las carpetas de animales")
    parser.add_argument(
        "--out-dir", default="results/Mapas_por_tipo",
        help="Carpeta base; se crean subcarpetas por animal y metadata de sesión",
    )
    parser.add_argument(
        "--exclude-day", action="append", default=[],
        help="Excluir una etiqueta sess.day_name exacta; se puede repetir (p. ej. HabL)",
    )
    args = parser.parse_args()
    process_all(
        data_dir=args.data_dir,
        out_dir=args.out_dir,
        animal=args.animal,
        session_file=args.session_file,
        exclude_day_names=args.exclude_day,
    )


if __name__ == "__main__":
    main()
