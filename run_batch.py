import argparse
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from core.data_loader import (
    get_common_neurons,
    list_sessions,
    load_session,
    mapping_columns,
)
from core.place_fields import compute_rate_map, draw_neuron_analysis, threshold_event_frames


def process_session(animal, session_file, data_dir=None, out_dir="results"):
    print(f"Procesando {animal} - {session_file}...", flush=True)
    data = load_session(animal, session_file, data_dir=data_dir)

    if data["type"] == "simple":
        subsessions = [data]
        common_neurons = np.arange(data["n_neurons"], dtype=int)
        columns_per_sub = [common_neurons.copy()]
        mapping = None
    else:
        subsessions = data["subsessions"]
        mapping = data["mapping"]
        common_neurons = get_common_neurons(mapping)
        columns_per_sub = [
            mapping_columns(mapping, common_neurons, index, sub["S"].shape[1])
            for index, sub in enumerate(subsessions)
        ]

    if len(common_neurons) == 0:
        print(f"  Sin neuronas comunes: {animal} - {session_file}", flush=True)
        return

    session_out = os.path.join(out_dir, animal, session_file.replace(".mat", ""))
    os.makedirs(session_out, exist_ok=True)
    common_index = {int(global_id): index for index, global_id in enumerate(common_neurons)}

    for global_id in common_neurons:
        fig, axes = plt.subplots(
            len(subsessions), 2, figsize=(10, 4 * len(subsessions))
        )
        axes = np.atleast_2d(axes)
        neuron_data = {}

        for sub_index, sub in enumerate(subsessions):
            column_index = common_index[int(global_id)]
            local_id = int(columns_per_sub[sub_index][column_index])
            if mapping is None:
                local_original_id = int(global_id)
            else:
                local_original_id = mapping[global_id, sub_index]

            signal = sub["S"][:, local_id]
            events, threshold = threshold_event_frames(signal)
            rate_map, x_bins, y_bins = compute_rate_map(
                sub["x"], sub["y"], sub["t"], events,
                bin_size=2.5,
                sigma=2.5,
            )
            extent = [x_bins[0], x_bins[-1], y_bins[0], y_bins[-1]]

            neuron_data[f"sub_{sub_index}_rate_map"] = rate_map
            neuron_data[f"sub_{sub_index}_events"] = events
            neuron_data[f"sub_{sub_index}_threshold_3sd"] = threshold

            ax_traj = axes[sub_index, 0]
            ax_map = axes[sub_index, 1]
            image = draw_neuron_analysis(
                ax_traj,
                ax_map,
                sub["x"],
                sub["y"],
                events,
                rate_map,
                extent,
                neuron_id=(
                    f"Global {global_id} "
                    f"(Local {local_original_id}->Col {local_id})"
                ),
            )
            ax_traj.set_title(f"Sub {sub_index + 1} | {ax_traj.get_title()}")
            ax_map.set_title(f"Sub {sub_index + 1} | {ax_map.get_title()}")
            fig.colorbar(image, ax=ax_map, fraction=0.046, pad=0.04).set_label(
                "eventos/s"
            )

        plt.tight_layout()
        png_path = os.path.join(session_out, f"neuron_{int(global_id):04d}.png")
        plt.savefig(png_path, dpi=100)
        plt.close(fig)

        npz_path = os.path.join(session_out, f"neuron_{int(global_id):04d}.npz")
        np.savez_compressed(npz_path, **neuron_data)

    print(
        f"Listo: {len(common_neurons)} neuronas en {session_out}",
        flush=True,
    )


def process_all(data_dir=None, out_dir="results", animal=None, session_file=None):
    sessions = [
        session for session in list_sessions(data_dir=data_dir)
        if (animal is None or session["animal"] == animal)
        and (session_file is None or session["session_file"] == session_file)
    ]
    if not sessions:
        raise SystemExit("No se encontraron sesiones que coincidan con los filtros.")

    for session in sessions:
        try:
            process_session(
                session["animal"],
                session["session_file"],
                data_dir=data_dir,
                out_dir=out_dir,
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
    parser.add_argument("--out-dir", default="results", help="Carpeta de salida")
    args = parser.parse_args()
    process_all(
        data_dir=args.data_dir,
        out_dir=args.out_dir,
        animal=args.animal,
        session_file=args.session_file,
    )


if __name__ == "__main__":
    main()
