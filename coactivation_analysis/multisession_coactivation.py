import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy

# Asegurar que el repo esté en el path
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from core.data_loader import (
    get_common_neurons,
    list_sessions,
    load_session,
    mapping_columns,
)
from core.place_fields import threshold_event_frames
from coactivation_analysis.metrics import (
    compute_pairwise_correlation,
    compute_jaccard_cooccurrence,
    find_coactive_assemblies
)


def get_subsession_labels(sess, n_subs):
    """Obtiene etiquetas descriptivas para las subsesiones."""
    default_labels = [
        "1. Open Field (AV1)",
        "2. 4 Objetos (Sample)",
        "3. 2 Objetos (Test)",
        "4. Open Field (AV2)"
    ]
    if sess is not None and hasattr(sess, 'names'):
        raw_names = sess.names
        if isinstance(raw_names, (list, np.ndarray)) and len(raw_names) == n_subs:
            labels = []
            for i, name in enumerate(raw_names):
                sname = str(name).strip(" []'")
                if not sname:
                    sname = default_labels[i] if i < len(default_labels) else f"Sub {i+1}"
                labels.append(f"Sub {i+1}: {sname}")
            return labels

    return default_labels[:n_subs] if n_subs <= len(default_labels) else [f"Sub {i+1}" for i in range(n_subs)]


def analyze_merged_coactivation(animal, session_file, out_base_dir=None, metric='pearson'):
    """
    Analiza la co-activación de neuronas comunes a lo largo de las 4 sub-sesiones
    de un archivo merged (evaluando la evolución de clústeres con objetos).
    """
    if out_base_dir is None:
        out_base_dir = os.path.join(REPO_ROOT, "coactivation_analysis", "results")

    print(f"\n========================================================")
    print(f"Análisis Multi-subsesión (Evolución de Clústeres con Objetos)")
    print(f"Animal: {animal} | Archivo: {session_file} | Métrica: {metric}")
    print(f"========================================================")

    data = load_session(animal, session_file)
    if data['type'] != 'merged':
        raise ValueError(f"La sesión {session_file} no es merged. Tipo: {data['type']}")

    mapping = data['mapping']
    subsessions = data['subsessions']
    n_subs = data['n_subsessions']
    sess = data['sess']

    common_neurons = get_common_neurons(mapping)
    n_common = len(common_neurons)
    print(f"Neuronas comunes presentes en todas las {n_subs} subsesiones: {n_common}")

    if n_common < 3:
        print("[AVISO] Muy pocas neuronas comunes (< 3). Omitiendo análisis multi-subsesión.")
        return None

    # Mapear columnas locales para cada subsesión
    columns_per_sub = [
        mapping_columns(mapping, common_neurons, i, subsessions[i]['S'].shape[1])
        for i in range(n_subs)
    ]

    sub_labels = get_subsession_labels(sess, n_subs)

    # Extraer la actividad alineada para las neuronas comunes
    aligned_activities = []
    corr_matrices = []

    for i in range(n_subs):
        local_cols = columns_per_sub[i]
        # Actividad S para neuronas comunes
        S_common = subsessions[i]['S'][:, local_cols]
        aligned_activities.append(S_common)

        if metric == 'pearson':
            R = compute_pairwise_correlation(S_common, method='pearson')
        elif metric == 'spearman':
            R = compute_pairwise_correlation(S_common, method='spearman')
        elif metric == 'jaccard':
            # Binarizar eventos
            ev_bin = np.zeros_like(S_common)
            for j in range(n_common):
                col = S_common[:, j]
                if np.std(col) > 0:
                    ev, _ = threshold_event_frames(col, threshold_sd=3.0)
                    ev_bin[:, j] = ev
            R = compute_jaccard_cooccurrence(ev_bin)
        else:
            raise ValueError(f"Métrica desconocida: {metric}")

        corr_matrices.append(R)

    corr_matrices = np.array(corr_matrices) # forma: (n_subs, n_common, n_common)

    # 1. Definir orden de clústeres basado en la Subsesión 1 (Línea de base / Open Field inicial)
    base_matrix = corr_matrices[0]
    base_order, base_labels, _ = find_coactive_assemblies(base_matrix, distance_threshold=0.5)

    # 2. Definir orden de clústeres basado en la correlación promedio
    mean_matrix = np.nanmean(corr_matrices, axis=0)
    mean_order, mean_labels, _ = find_coactive_assemblies(mean_matrix, distance_threshold=0.5)

    # Crear directorio de salida
    sess_id = f"{animal}_{os.path.splitext(session_file)[0]}"
    out_dir = os.path.join(out_base_dir, sess_id, "multisession_evolution")
    os.makedirs(out_dir, exist_ok=True)

    # 3. Guardar matrices ordenadas
    np.save(os.path.join(out_dir, "corr_matrices_raw.npy"), corr_matrices)
    np.save(os.path.join(out_dir, "common_neurons_global_ids.npy"), common_neurons)
    np.save(os.path.join(out_dir, "base_order.npy"), base_order)

    # 4. Calcular estabilidad / similitud entre las matrices de correlación (matriz 4x4)
    triu_idx = np.triu_indices(n_common, k=1)
    matrix_sim = np.ones((n_subs, n_subs), dtype=float)
    for i in range(n_subs):
        for j in range(i + 1, n_subs):
            vals_i = corr_matrices[i][triu_idx]
            vals_j = corr_matrices[j][triu_idx]
            valid = np.isfinite(vals_i) & np.isfinite(vals_j)
            if np.sum(valid) > 2:
                r_sim = float(np.corrcoef(vals_i[valid], vals_j[valid])[0, 1])
            else:
                r_sim = np.nan
            matrix_sim[i, j] = r_sim
            matrix_sim[j, i] = r_sim

    sim_df = pd.DataFrame(matrix_sim, index=[f"Sub {i+1}" for i in range(n_subs)],
                          columns=[f"Sub {i+1}" for i in range(n_subs)])
    sim_df.to_csv(os.path.join(out_dir, "matrix_similarity_stability.csv"))
    print("\n--- Matriz de Similitud entre Subsesiones (Estabilidad de Co-activación) ---")
    print(sim_df.round(3))

    # 5. Generar Figura 1: Evolución de las Matrices con Orden FIJO basado en Subsesión 1 (Baseline)
    fig, axes = plt.subplots(1, n_subs, figsize=(5 * n_subs, 4.5), dpi=150)
    if n_subs == 1:
        axes = [axes]

    vmin = 0.0 if metric == 'jaccard' else -1.0
    vmax = 1.0 if metric == 'jaccard' else 1.0
    cmap = 'viridis' if metric == 'jaccard' else 'coolwarm'

    for i in range(n_subs):
        ax = axes[i]
        mat = corr_matrices[i][np.ix_(base_order, base_order)]
        im = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax, origin='upper', aspect='equal')
        ax.set_title(f"{sub_labels[i]}", fontsize=11, fontweight='bold', pad=8)
        ax.set_xlabel("Neuronas (Orden Sub 1)")
        if i == 0:
            ax.set_ylabel("Neuronas (Orden Sub 1)")
        else:
            ax.set_yticks([])

    # Colorbar compartida
    cbar = fig.colorbar(im, ax=axes, shrink=0.75, pad=0.02)
    cbar.set_label(f"Co-activación ({metric.capitalize()})", fontsize=11)

    fig.suptitle(
        f"Evolución de Clusters de Co-activación (Mismo orden fijado en Subsesión 1 - Baseline)\n"
        f"{animal} - {session_file} ({n_common} neuronas comunes)",
        fontsize=13, y=1.03
    )

    plot_fixed_path = os.path.join(out_dir, f"coactivation_evolution_fixed_order_{metric}.png")
    plt.savefig(plot_fixed_path, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"\n[OK] Gráfico de evolución fijada guardado en:\n  {plot_fixed_path}")

    # 6. Generar Figura 2: Re-clustering Independiente por Subsesión
    # (Para ver si en la presencia de objetos se forman ensambles totalmente nuevos)
    fig2, axes2 = plt.subplots(1, n_subs, figsize=(5 * n_subs, 4.5), dpi=150)
    if n_subs == 1:
        axes2 = [axes2]

    for i in range(n_subs):
        ax = axes2[i]
        ord_i, _, _ = find_coactive_assemblies(corr_matrices[i], distance_threshold=0.5)
        mat_i = corr_matrices[i][np.ix_(ord_i, ord_i)]
        im2 = ax.imshow(mat_i, cmap=cmap, vmin=vmin, vmax=vmax, origin='upper', aspect='equal')
        ax.set_title(f"{sub_labels[i]}\n(Re-clusterizado)", fontsize=10, pad=8)
        ax.set_xlabel("Neuronas (Orden Propio)")
        if i == 0:
            ax.set_ylabel("Neuronas")
        else:
            ax.set_yticks([])

    cbar2 = fig2.colorbar(im2, ax=axes2, shrink=0.75, pad=0.02)
    cbar2.set_label(f"Co-activación ({metric.capitalize()})", fontsize=11)

    fig2.suptitle(
        f"Estructura Modular Independiente por Condición (Objetos vs Open Field)\n"
        f"{animal} - {session_file} ({n_common} neuronas comunes)",
        fontsize=13, y=1.03
    )

    plot_recluster_path = os.path.join(out_dir, f"coactivation_recluster_per_session_{metric}.png")
    plt.savefig(plot_recluster_path, bbox_inches='tight', dpi=300)
    plt.close(fig2)
    print(f"[OK] Gráfico de re-clustering independiente guardado en:\n  {plot_recluster_path}")

    # 7. Generar Figura 3: Matriz de Similitud 4x4
    fig3, ax3 = plt.subplots(figsize=(5, 4.2), dpi=150)
    im3 = ax3.imshow(matrix_sim, cmap='viridis', vmin=0, vmax=1.0)
    ax3.set_xticks(range(n_subs))
    ax3.set_yticks(range(n_subs))
    ax3.set_xticklabels([f"Sub {i+1}" for i in range(n_subs)], rotation=45, ha='right')
    ax3.set_yticklabels([f"Sub {i+1}" for i in range(n_subs)])
    for r in range(n_subs):
        for c in range(n_subs):
            val = matrix_sim[r, c]
            ax3.text(c, r, f"{val:.2f}", ha='center', va='center',
                     color='white' if val < 0.6 else 'black', fontweight='bold')
    ax3.set_title("Estabilidad de Co-activación\n(Correlación entre matrices)", fontsize=11)
    fig3.colorbar(im3, ax=ax3, shrink=0.8)
    plt.tight_layout()
    plot_sim_path = os.path.join(out_dir, f"matrix_similarity_heatmap.png")
    plt.savefig(plot_sim_path, bbox_inches='tight', dpi=300)
    plt.close(fig3)
    print(f"[OK] Heatmap de estabilidad guardado en:\n  {plot_sim_path}")

    return {
        'n_common': n_common,
        'common_neurons': common_neurons,
        'corr_matrices': corr_matrices,
        'base_order': base_order,
        'matrix_similarity': matrix_sim,
        'sub_labels': sub_labels,
        'out_dir': out_dir,
        'plot_fixed_path': plot_fixed_path,
        'plot_recluster_path': plot_recluster_path,
        'plot_sim_path': plot_sim_path
    }


def main():
    parser = argparse.ArgumentParser(
        description="Evolución de Matrices de Co-activación en Sesiones Merged"
    )
    parser.add_argument("--animal", type=str, help="Animal ID (ej. R004)")
    parser.add_argument("--session", type=str, help="Archivo de sesión merged")
    parser.add_argument(
        "--metric",
        type=str,
        default="pearson",
        choices=["pearson", "spearman", "jaccard"],
    )
    args = parser.parse_args()

    sessions = [
        session for session in list_sessions()
        if session["session_file"].endswith("_merged.mat")
        and (args.animal is None or session["animal"] == args.animal)
        and (args.session is None or session["session_file"] == args.session)
    ]
    if not sessions:
        raise SystemExit("No se encontraron archivos merged que coincidan con los filtros.")

    processed = 0
    for session in sessions:
        try:
            data = load_session(session["animal"], session["session_file"])
        except Exception as error:
            print(
                f"[OMITIDA] {session['animal']} / {session['session_file']}: {error}"
            )
            continue
        if data["type"] != "merged":
            print(
                f"[OMITIDA] {session['animal']} / {session['session_file']}: "
                "no es una sesión merged válida"
            )
            continue
        analyze_merged_coactivation(
            session["animal"], session["session_file"], metric=args.metric
        )
        processed += 1

    if processed == 0:
        raise SystemExit("No hubo sesiones merged válidas para procesar.")


if __name__ == "__main__":
    main()
