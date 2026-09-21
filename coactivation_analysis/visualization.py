import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_coactivation_matrix(corr_matrix, order=None, labels=None, title="Matriz de Co-activación Neuronal", save_path=None, vmin=None, vmax=None):
    """
    Grafica la matriz de correlación / co-activación par a par.

    Parámetros
    ----------
    corr_matrix : ndarray, forma (n_neurons, n_neurons)
    order : ndarray, opcional
        Índices de reordenamiento (p. ej. por clustering jerárquico).
    labels : ndarray, opcional
        Etiquetas de ensambles / clústeres.
    title : str
    save_path : str, opcional
    """
    matrix_to_plot = corr_matrix.copy()
    if order is not None:
        matrix_to_plot = matrix_to_plot[np.ix_(order, order)]

    fig, ax = plt.subplots(figsize=(8, 7), dpi=150)

    cmap = 'viridis' if np.nanmin(corr_matrix) >= 0 else 'coolwarm'
    if vmin is None:
        vmin = 0.0 if np.nanmin(corr_matrix) >= 0 else -1.0
    if vmax is None:
        vmax = 1.0

    im = ax.imshow(matrix_to_plot, cmap=cmap, origin='upper', aspect='equal', vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Correlación / Co-activación")

    ax.set_title(title, fontsize=14, pad=12)
    ax.set_xlabel("Neuronas" + (" (Reordenadas)" if order is not None else ""))
    ax.set_ylabel("Neuronas" + (" (Reordenadas)" if order is not None else ""))

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        plt.close(fig)
    return fig, ax


def plot_coactive_pairs_raster(activity, top_pairs, t=None, title="Rasters de Actividad de Pares Co-activos", save_path=None, max_frames=1000):
    """
    Muestra trazas de tiempo alineadas para los pares de neuronas con mayor co-disparo.

    Parámetros
    ----------
    activity : ndarray, forma (n_frames, n_neurons)
    top_pairs : list of tuples [(n1, n2, score), ...]
    t : ndarray, opcional
        Vector de tiempo en segundos.
    title : str
    save_path : str, opcional
    max_frames : int
        Límite de frames para visualización clara.
    """
    n_frames, n_neurons = activity.shape
    display_frames = min(n_frames, max_frames)

    if t is None:
        t = np.arange(display_frames) * 0.05
    else:
        t = t[:display_frames]

    n_pairs = len(top_pairs)
    if n_pairs == 0:
        return None, None

    fig, axes = plt.subplots(n_pairs, 1, figsize=(10, 2.5 * n_pairs), sharex=True, dpi=150)
    if n_pairs == 1:
        axes = [axes]

    for idx, (n1, n2, score) in enumerate(top_pairs):
        ax = axes[idx]
        sig1 = activity[:display_frames, n1]
        sig2 = activity[:display_frames, n2]

        ax.plot(t, sig1, label=f"Neurona {n1}", color='#1f77b4', alpha=0.85, linewidth=1.2)
        ax.plot(t, sig2, label=f"Neurona {n2}", color='#ff7f0e', alpha=0.85, linewidth=1.2)

        ax.set_title(f"Par #{idx+1}: Neurona {n1} y Neurona {n2} (Score: {score:.3f})", fontsize=11)
        ax.set_ylabel("Actividad")
        ax.legend(loc='upper right', fontsize=9)
        ax.grid(True, linestyle='--', alpha=0.4)

    axes[-1].set_xlabel("Tiempo (s)")
    fig.suptitle(title, fontsize=14, y=1.02)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        plt.close(fig)
    return fig, axes


def plot_coactivity_distribution(observed_corr, shuffled_mean=None, z_scores=None, title="Distribución de Correlaciones de Co-activación", save_path=None):
    """
    Grafica la distribución de correlaciones observadas par a par y las compara con el modelo nulo.
    """
    # Extraer triángulo superior (excluyendo diagonal)
    n_neurons = observed_corr.shape[0]
    triu_idx = np.triu_indices(n_neurons, k=1)

    obs_vals = observed_corr[triu_idx]
    obs_vals = obs_vals[np.isfinite(obs_vals)]

    fig, axes = plt.subplots(1, 2 if z_scores is not None else 1, figsize=(12 if z_scores is not None else 6, 4.5), dpi=150)
    if z_scores is None:
        ax1 = axes
    else:
        ax1, ax2 = axes

    sns.histplot(obs_vals, ax=ax1, kde=True, color='teal', stat='density', bins=30, label='Observada')

    if shuffled_mean is not None:
        shuff_vals = shuffled_mean[triu_idx]
        shuff_vals = shuff_vals[np.isfinite(shuff_vals)]
        sns.histplot(shuff_vals, ax=ax1, kde=True, color='gray', stat='density', bins=30, label='Modelo Nulo (Shift)', alpha=0.5)

    ax1.set_title("Correlación Par a Par", fontsize=12)
    ax1.set_xlabel("Correlación ($r$ / Jaccard)")
    ax1.set_ylabel("Densidad")
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.3)

    if z_scores is not None:
        z_vals = z_scores[triu_idx]
        z_vals = z_vals[np.isfinite(z_vals)]
        sns.histplot(z_vals, ax=ax2, kde=True, color='darkred', stat='density', bins=30)
        ax2.axvline(2.0, color='red', linestyle='--', label='Significativo (z = +2)')
        ax2.set_title("Distribución de Z-scores", fontsize=12)
        ax2.set_xlabel("Z-score vs. Modelo Nulo")
        ax2.set_ylabel("Densidad")
        ax2.legend()
        ax2.grid(True, linestyle='--', alpha=0.3)

    plt.suptitle(title, fontsize=14, y=1.03)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        plt.close(fig)
    return fig, axes


def plot_cross_correlation_lags(cross_corr, lags, top_pairs, frame_dt=0.05, title="Correlación Cruzada con Lags Temporales", save_path=None):
    """
    Muestra las curvas de correlación cruzada en función del lag de tiempo (segundos) para los pares destacados.
    """
    n_pairs = len(top_pairs)
    if n_pairs == 0:
        return None, None

    lag_times_s = lags * frame_dt

    fig, ax = plt.subplots(figsize=(8, 5), dpi=150)

    for idx, (n1, n2, score) in enumerate(top_pairs):
        cc_curve = cross_corr[n1, n2, :]
        ax.plot(lag_times_s, cc_curve, label=f"Par {n1}-{n2} (r={score:.2f})", linewidth=1.8)

    ax.axvline(0.0, color='black', linestyle=':', alpha=0.7, label='Lag cero')
    ax.set_title(title, fontsize=13)
    ax.set_xlabel("Lag de tiempo (segundos)")
    ax.set_ylabel("Correlación cruzada")
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, linestyle='--', alpha=0.4)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        plt.close(fig)
    return fig, ax
