"""generar_diagnostico_dlc_excluidos.py

Genera una imagen diagnóstica con la trayectoria del animal (sin actividad neuronal)
para las sesiones con problemas conocidos en la inferencia de DeepLabCut:
- R004 - T3_XsS_CNO (2026_06_24)
- R005 - T3_XsS_VEH (2026_06_24)

Estas sesiones presentan saltos bruscos y repetitivos de posicion hacia un punto
determinado del campo o fuera de los limites de la arena. Para evitar contaminar
el batch de place fields con mapas anomalos, se excluyen del batch principal
(run_batch.py) y se documentan mediante este script, conservando exactamente la
misma jerarquia de carpetas.
"""

import argparse
import os
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io

from core.data_loader import get_data_dir, list_sessions, load_session


ARENA_RADIUS_CM = 25.0


def _session_day_name(animal, session_file, data_dir):
    path = Path(data_dir) / animal / session_file
    metadata = scipy.io.loadmat(
        str(path), variable_names=["sess"], squeeze_me=True,
        struct_as_record=False,
    )
    sess = metadata.get("sess")
    return str(getattr(sess, "day_name", "")) if sess is not None else ""


def _session_category(day_name):
    """Estructura de carpetas derivada de sess.day_name identica a run_batch.py."""
    if day_name == "HabL":
        return ("HabL",)
    if re.fullmatch(r"HabC\d*", day_name):
        return (day_name,)
    match = re.fullmatch(r"T\d+_(SD|XsS)_(CNO|VEH)_?", day_name)
    if match:
        task, treatment = match.groups()
        return ("Objetos", task, treatment)
    return ("Otros",)


def _phase_names(day_name, n_subsessions):
    if n_subsessions == 1:
        return [day_name or "Registro"]
    if day_name == "HabL":
        return ["OF1", "OF2", "OF3", "OF4"]
    if "_SD_" in day_name or "_XsS_" in day_name:
        return ["OF1", "SAMPLE", "TEST", "OF2"]
    return [f"Sub {index + 1}" for index in range(n_subsessions)]


def is_dlc_defect_target(animal, session_file, day_name=""):
    """Detecta exclusivamente R004 XsS CNO y R005 XsS VEH."""
    a = str(animal).upper()
    s = str(session_file).upper()
    d = str(day_name).upper()

    is_r004 = "R004" in a or "M4" in a
    is_r004_target = is_r004 and (
        ("XSS" in d and "CNO" in d) or
        ("XSS" in s and "CNO" in s) or
        ("2026_06_24" in s and ("CNO" in d or "CNO" in s or not d))
    )

    is_r005 = "R005" in a or "M5" in a
    is_r005_target = is_r005 and (
        ("XSS" in d and "VEH" in d) or
        ("XSS" in s and "VEH" in s) or
        ("2026_06_24" in s and ("VEH" in d or "VEH" in s or not d))
    )

    return is_r004_target or is_r005_target


def render_trajectory_diagnostic(animal, session_file, day_name, subsessions, session_out):
    """Dibuja un panel 1x4 con la trayectoria del animal y resalta anomalias de tracking."""
    phases = _phase_names(day_name, len(subsessions))
    n_phases = len(subsessions)

    fig, axes = plt.subplots(1, n_phases, figsize=(3.8 * n_phases, 4.4), squeeze=False)
    theta = np.linspace(0, 2 * np.pi, 200)
    arena_x = ARENA_RADIUS_CM * np.cos(theta)
    arena_y = ARENA_RADIUS_CM * np.sin(theta)

    for i, sub in enumerate(subsessions):
        ax = axes[0, i]
        ax.set_aspect("equal")
        ax.set_xlim(-ARENA_RADIUS_CM - 6, ARENA_RADIUS_CM + 6)
        ax.set_ylim(-ARENA_RADIUS_CM - 6, ARENA_RADIUS_CM + 6)
        ax.set_xticks((-25, 0, 25))
        ax.set_yticks((-25, 0, 25))
        ax.tick_params(labelsize=7)

        # Limite arena circular
        ax.plot(arena_x, arena_y, color="black", linestyle="--", linewidth=1.0, alpha=0.7, zorder=1)

        x = np.asarray(sub["x"], dtype=float)
        y = np.asarray(sub["y"], dtype=float)
        valid = np.isfinite(x) & np.isfinite(y)
        x_val, y_val = x[valid], y[valid]

        r = np.sqrt(x_val**2 + y_val**2)
        out_arena = r > (ARENA_RADIUS_CM + 1.0)
        n_out = int(np.sum(out_arena))
        pct_out = (n_out / len(x_val) * 100) if len(x_val) > 0 else 0.0

        # Traza trayectoria
        ax.plot(x_val, y_val, color="0.45", linewidth=0.4, alpha=0.75, rasterized=True, zorder=2)

        # Resaltar puntos anomalos o fuera de arena
        if n_out > 0:
            ax.scatter(
                x_val[out_arena], y_val[out_arena],
                color="#d62728", s=6, alpha=0.9,
                label=f"Fuera arena ({n_out})", zorder=3, rasterized=True
            )

        ax.set_title(
            f"{phases[i]}\n{len(x_val)} frames | {n_out} ({pct_out:.1f}%) fuera arena",
            fontsize=9,
        )

    fig.suptitle(
        f"{animal} | {day_name} ({session_file}) - SESIÓN EXCLUIDA DEL BATCH DE MAPAS\n"
        f"Motivo: Artefactos severos de tracking en DeepLabCut (saltos bruscos hacia punto fijo / fuera de límites)",
        fontsize=10, fontweight="bold", color="#b22222", y=0.98,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))

    img_path = os.path.join(session_out, "trayectoria_anomala_deeplabcut.png")
    fig.savefig(img_path, dpi=140, facecolor="white")
    plt.close(fig)
    print(f"  [OK] Guardada trayectoria diagnostica en: {img_path}", flush=True)

    # Dejar archivo de texto explicativo en la carpeta
    readme_path = os.path.join(session_out, "LEEME_EXCLUSION_DLC.txt")
    with open(readme_path, "w", encoding="utf-8") as f:
        f.write(
            f"SESION EXCLUIDA DE LA GENERACION DE MAPAS ESPACIALES (PLACE FIELDS)\n"
            f"{'='*70}\n\n"
            f"Animal: {animal}\n"
            f"Archivo: {session_file}\n"
            f"Day Name: {day_name}\n\n"
            f"MOTIVO DE LA EXCLUSION:\n"
            f"Esta sesion presenta artefactos graves en la inferencia de tracking de DeepLabCut (DLC),\n"
            f"evidenciados por saltos bruscos y anomalos de posicion hacia un punto determinado del campo\n"
            f"o fuera del perimetro circular de la arena.\n\n"
            f"Para evitar contaminar el analisis con place fields espurios causados por errores del modelo\n"
            f"de tracking conductual, esta sesion fue omitida del batch principal de generacion de figuras.\n\n"
            f"DIAGNOSTICO VISUAL:\n"
            f"Consulte el archivo 'trayectoria_anomala_deeplabcut.png' en esta misma carpeta para\n"
            f"inspeccionar la trayectoria del animal y los saltos anomalos (sin actividad neuronal).\n"
        )
    print(f"  [OK] Guardado LEEME explicativo en: {readme_path}", flush=True)


def process_excluded_sessions(data_dir=None, out_dir="results/Mapas_por_tipo"):
    data_dir = data_dir or get_data_dir()
    sessions = list_sessions(data_dir=data_dir)
    target_sessions = []

    for s in sessions:
        try:
            day_name = _session_day_name(s["animal"], s["session_file"], data_dir)
        except Exception:
            day_name = ""
        if is_dlc_defect_target(s["animal"], s["session_file"], day_name):
            target_sessions.append((s, day_name))

    if not target_sessions:
        print(f"No se encontraron sesiones objetivo (R004 XsS CNO o R005 XsS VEH) en {data_dir}")
        return

    print(f"Procesando {len(target_sessions)} sesiones con anomalias de DLC...")
    for s, day_name in target_sessions:
        animal = s["animal"]
        session_file = s["session_file"]
        category_parts = _session_category(day_name)
        session_out = os.path.join(
            out_dir, animal, *category_parts,
            session_file.replace(".mat", ""),
        )
        os.makedirs(session_out, exist_ok=True)
        print(f"Procesando {animal} - {session_file} ({day_name}) -> {session_out}...")
        data = load_session(animal, session_file, data_dir=data_dir)
        subsessions = [data] if data["type"] == "simple" else data["subsessions"]
        render_trajectory_diagnostic(animal, session_file, day_name, subsessions, session_out)


def main():
    parser = argparse.ArgumentParser(
        description="Genera la figura de diagnostico de trayectoria para sesiones excluidas por artefactos de DLC."
    )
    parser.add_argument("--data-dir", default=None, help="Carpeta raiz de datos (por defecto MINISCOPE_DATA_DIR)")
    parser.add_argument(
        "--out-dir", default="results/Mapas_por_tipo",
        help="Carpeta base de salida (identica a la de run_batch.py)",
    )
    args = parser.parse_args()
    process_excluded_sessions(data_dir=args.data_dir, out_dir=args.out_dir)


if __name__ == "__main__":
    main()
