import os
import sys
import argparse
import numpy as np
import pandas as pd
import scipy.io

# Asegurar que el directorio raíz del repositorio esté en sys.path
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from core.data_loader import get_data_dir, list_sessions, load_session
from core.place_fields import threshold_event_frames
from coactivation_analysis.metrics import (
    compute_pairwise_correlation,
    compute_jaccard_cooccurrence,
    compute_cross_correlation_lags,
    compute_shuffled_significance,
    find_coactive_assemblies
)
from coactivation_analysis.visualization import (
    plot_coactivation_matrix,
    plot_coactive_pairs_raster,
    plot_coactivity_distribution,
    plot_cross_correlation_lags
)


def load_continuous_C(animal, session_file):
    """Auxiliar para cargar la señal de fluorescencia continua C si existe."""
    fpath = os.path.join(get_data_dir(), animal, session_file)
    try:
        d = scipy.io.loadmat(fpath, squeeze_me=True, struct_as_record=False)
        if 'act' in d and hasattr(d['act'], 'C'):
            return np.asarray(d['act'].C, dtype=float)
        elif 'C' in d:
            return np.asarray(d['C'], dtype=float)
    except Exception:
        pass
    return None


def get_top_pairs(corr_matrix, top_k=5):
    """Devuelve los top K pares de neuronas distintas con mayor correlación/co-activación."""
    n_neurons = corr_matrix.shape[0]
    if n_neurons < 2:
        return []

    triu_i, triu_j = np.triu_indices(n_neurons, k=1)
    values = corr_matrix[triu_i, triu_j]

    # Filtrar NaN / Infs
    valid = np.isfinite(values)
    triu_i, triu_j, values = triu_i[valid], triu_j[valid], values[valid]

    if len(values) == 0:
        return []

    top_indices = np.argsort(values)[::-1][:top_k]
    return [(int(triu_i[idx]), int(triu_j[idx]), float(values[idx])) for idx in top_indices]


def analyze_single_activity(activity, t, out_dir, session_name, signal_label="S_deconvolved", n_shuffles=50):
    """
    Ejecuta el pipeline completo de co-activación para una matriz de actividad (frames x neuronas).
    """
    os.makedirs(out_dir, exist_ok=True)
    n_frames, n_neurons = activity.shape
    print(f"  -> Procesando {session_name} [{signal_label}] ({n_frames} frames, {n_neurons} neuronas)...")

    if n_neurons < 2:
        print(f"  [AVISO] Insuficientes neuronas ({n_neurons}) en {session_name}. Omitiendo.")
        return

    # 1. Binarización de eventos (S > 3*SD)
    events_bin = np.zeros_like(activity)
    for j in range(n_neurons):
        col = activity[:, j]
        if np.std(col) > 0:
            ev, _ = threshold_event_frames(col, threshold_sd=3.0)
            events_bin[:, j] = ev

    # 2. Matrices de correlación
    pearson_corr = compute_pairwise_correlation(activity, method='pearson')
    spearman_corr = compute_pairwise_correlation(activity, method='spearman')
    jaccard_matrix = compute_jaccard_cooccurrence(events_bin)

    # 3. Test de significancia con modelo nulo (Shuffling)
    shuff_results = compute_shuffled_significance(
        activity, n_shuffles=n_shuffles, min_shift_s=30.0, frame_dt=0.05, method='pearson'
    )

    # 4. Correlación cruzada con lags temporales (+-20 frames = +-1s)
    cross_corr, lags = compute_cross_correlation_lags(activity, max_lag_frames=20)

    # 5. Clustering jerárquico de ensambles co-activos
    order, labels, _ = find_coactive_assemblies(pearson_corr, distance_threshold=0.5)

    # 6. Top pares co-activos
    top_pearson = get_top_pairs(pearson_corr, top_k=5)
    top_jaccard = get_top_pairs(jaccard_matrix, top_k=5)

    # --- Guardar Resultados Cuantitativos ---
    np.save(os.path.join(out_dir, f"{session_name}_{signal_label}_pearson.npy"), pearson_corr)
    np.save(os.path.join(out_dir, f"{session_name}_{signal_label}_jaccard.npy"), jaccard_matrix)
    np.save(os.path.join(out_dir, f"{session_name}_{signal_label}_zscores.npy"), shuff_results['z_scores'])

    # Guardar CSV de pares destacados
    top_df = pd.DataFrame([
        {
            'rank': idx + 1,
            'neuron_1': n1,
            'neuron_2': n2,
            'pearson_r': score,
            'jaccard_index': jaccard_matrix[n1, n2],
            'z_score': shuff_results['z_scores'][n1, n2],
            'p_value': shuff_results['p_values'][n1, n2]
        }
        for idx, (n1, n2, score) in enumerate(top_pearson)
    ])
    top_df.to_csv(os.path.join(out_dir, f"{session_name}_{signal_label}_top_coactive_pairs.csv"), index=False)

    # --- Guardar Figuras / Gráficos ---
    plot_coactivation_matrix(
        pearson_corr, order=order,
        title=f"Matriz Co-activación Pearson ({session_name} - {signal_label})",
        save_path=os.path.join(out_dir, f"{session_name}_{signal_label}_matrix_pearson.png")
    )

    plot_coactivation_matrix(
        jaccard_matrix, order=order, vmin=0, vmax=float(np.nanmax(jaccard_matrix) or 1.0),
        title=f"Matriz Co-ocurrencia Jaccard ({session_name} - {signal_label})",
        save_path=os.path.join(out_dir, f"{session_name}_{signal_label}_matrix_jaccard.png")
    )

    plot_coactivity_distribution(
        pearson_corr, shuffled_mean=shuff_results['shuffled_mean'], z_scores=shuff_results['z_scores'],
        title=f"Distribución Correlaciones vs. Nulo ({session_name})",
        save_path=os.path.join(out_dir, f"{session_name}_{signal_label}_distribution.png")
    )

    if len(top_pearson) > 0:
        plot_coactive_pairs_raster(
            activity, top_pearson, t=t,
            title=f"Trazas Temporales - Top Pares Co-activos ({session_name})",
            save_path=os.path.join(out_dir, f"{session_name}_{signal_label}_rasters.png")
        )

        plot_cross_correlation_lags(
            cross_corr, lags, top_pearson, frame_dt=0.05,
            title=f"Correlación Cruzada (Lags) - Top Pares ({session_name})",
            save_path=os.path.join(out_dir, f"{session_name}_{signal_label}_crosscorr.png")
        )

    print(f"  [OK] Resultados guardados en: {out_dir}")


def process_session(session_info, out_base_dir, n_shuffles=50):
    animal = session_info['animal']
    session_file = session_info['session_file']
    sess_name = f"{animal}_{os.path.splitext(session_file)[0]}"

    print(f"\n==========================================")
    print(f"Cargando sesión: {animal} / {session_file}")
    print(f"==========================================")

    try:
        data = load_session(animal, session_file)
    except Exception as e:
        print(f"Error al cargar la sesión {session_file}: {e}")
        return

    sess_out_dir = os.path.join(out_base_dir, sess_name)

    if data['type'] == 'simple':
        S = data['S']
        t = data['t']
        analyze_single_activity(S, t, sess_out_dir, sess_name, signal_label="S_deconvolved", n_shuffles=n_shuffles)

    elif data['type'] == 'merged':
        n_subs = data['n_subsessions']
        sub_names = ['A', 'B', 'C', 'D']
        print(f"  Sesión merged con {n_subs} subsesiones.")

        for i, sub in enumerate(data['subsessions']):
            sub_id = sub_names[i] if i < len(sub_names) else f"sub{i+1}"
            sub_name = f"{sess_name}_{sub_id}"
            sub_out_dir = os.path.join(sess_out_dir, f"subsession_{sub_id}")

            S = sub['S']
            t = sub['t']
            analyze_single_activity(S, t, sub_out_dir, sub_name, signal_label="S_deconvolved", n_shuffles=n_shuffles)


def main():
    parser = argparse.ArgumentParser(description="Análisis de Co-activación y Correlación Neuronal")
    parser.add_argument("--animal", type=str, help="Filtro por animal (ej. R004)")
    parser.add_argument("--session", type=str, help="Filtro por nombre de archivo (ej. 2026_06_17_merged.mat)")
    parser.add_argument("--shuffles", type=int, default=50, help="Número de permutaciones para modelo nulo")
    parser.add_argument("--output_dir", type=str, default=os.path.join(REPO_ROOT, "coactivation_analysis", "results"))

    args = parser.parse_args()

    sessions = list_sessions()
    if not sessions:
        print("No se encontraron sesiones en la carpeta data/")
        return

    filtered_sessions = []
    for s in sessions:
        if args.animal and s['animal'] != args.animal:
            continue
        if args.session and s['session_file'] != args.session:
            continue
        filtered_sessions.append(s)

    print(f"Se procesarán {len(filtered_sessions)} sesión(es)...")

    for s_info in filtered_sessions:
        process_session(s_info, args.output_dir, n_shuffles=args.shuffles)

    print("\n¡Análisis de co-activación completado exitosamente!")


if __name__ == "__main__":
    main()
