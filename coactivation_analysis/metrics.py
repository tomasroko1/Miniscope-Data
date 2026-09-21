import numpy as np
from scipy import stats
from scipy.cluster import hierarchy

def compute_pairwise_correlation(activity, method='pearson'):
    """
    Calcula la matriz de correlación par a par entre neuronas.

    Parámetros
    ----------
    activity : ndarray, forma (n_frames, n_neurons)
        Matriz de actividad neuronal.
    method : str, opcional ('pearson', 'spearman', 'cosine')
        Método de correlación/similitud.

    Devuelve
    --------
    corr_matrix : ndarray, forma (n_neurons, n_neurons)
        Matriz simétrica de correlación.
    """
    activity = np.asarray(activity, dtype=float)
    if activity.ndim == 1:
        activity = activity[:, None]

    n_frames, n_neurons = activity.shape
    if n_neurons == 0:
        return np.empty((0, 0))

    if method == 'pearson':
        # Manejo eficiente de nan / varianza cero
        mean = np.nanmean(activity, axis=0, keepdims=True)
        std = np.nanstd(activity, axis=0, keepdims=True)

        norm_activity = np.where(std > 0, (activity - mean) / std, 0.0)
        # NaN a 0
        norm_activity = np.nan_to_num(norm_activity, nan=0.0)

        corr_matrix = (norm_activity.T @ norm_activity) / max(1, n_frames - 1)
        # Ajuste diagonal a 1 si hay varianza
        np.fill_diagonal(corr_matrix, 1.0)
        return np.clip(corr_matrix, -1.0, 1.0)

    elif method == 'spearman':
        # Rangos por columna
        ranks = np.zeros_like(activity)
        for i in range(n_neurons):
            col = activity[:, i]
            valid_mask = np.isfinite(col)
            if np.any(valid_mask):
                ranks[valid_mask, i] = stats.rankdata(col[valid_mask])
        return compute_pairwise_correlation(ranks, method='pearson')

    elif method == 'cosine':
        clean_act = np.nan_to_num(activity, nan=0.0)
        norms = np.linalg.norm(clean_act, axis=0, keepdims=True)
        norms = np.where(norms > 0, norms, 1.0)
        normalized = clean_act / norms
        cosine_sim = normalized.T @ normalized
        np.fill_diagonal(cosine_sim, 1.0)
        return np.clip(cosine_sim, 0.0, 1.0)

    else:
        raise ValueError(f"Método no soportado: {method}. Use 'pearson', 'spearman' o 'cosine'.")


def compute_jaccard_cooccurrence(events_bin):
    """
    Calcula el índice de Jaccard par a par para frames-evento binarios.

    J(i, j) = (Frames donde ambas están activas) / (Frames donde al menos una está activa)

    Parámetros
    ----------
    events_bin : ndarray, forma (n_frames, n_neurons)
        Matriz binaria (0 ó 1).

    Devuelve
    --------
    jaccard_matrix : ndarray, forma (n_neurons, n_neurons)
    """
    bin_act = (np.asarray(events_bin, dtype=float) > 0).astype(float)
    bin_act = np.nan_to_num(bin_act, nan=0.0)

    intersection = bin_act.T @ bin_act
    sum_active = np.sum(bin_act, axis=0, keepdims=True)
    union = sum_active.T + sum_active - intersection

    jaccard_matrix = np.divide(intersection, union, out=np.zeros_like(intersection), where=union > 0)
    np.fill_diagonal(jaccard_matrix, 1.0)
    return np.clip(jaccard_matrix, 0.0, 1.0)


def compute_cross_correlation_lags(activity, max_lag_frames=20):
    """
    Calcula la matriz de correlación cruzada para múltiples lags temporales.

    Parámetros
    ----------
    activity : ndarray, forma (n_frames, n_neurons)
    max_lag_frames : int
        Número máximo de frames de desplazamiento a izquierda y derecha.

    Devuelve
    --------
    cross_corr : ndarray, forma (n_neurons, n_neurons, 2 * max_lag + 1)
    lags : ndarray, forma (2 * max_lag + 1,)
    """
    activity = np.asarray(activity, dtype=float)
    activity = np.nan_to_num(activity, nan=0.0)

    # Normalizar por celula (media 0, std 1)
    mean = np.mean(activity, axis=0, keepdims=True)
    std = np.std(activity, axis=0, keepdims=True)
    std_safe = np.where(std > 0, std, 1.0)
    norm_act = (activity - mean) / std_safe

    n_frames, n_neurons = norm_act.shape
    lags = np.arange(-max_lag_frames, max_lag_frames + 1)
    n_lags = len(lags)

    cross_corr = np.zeros((n_neurons, n_neurons, n_lags), dtype=float)

    for idx, lag in enumerate(lags):
        if lag < 0:
            a = norm_act[-lag:, :]
            b = norm_act[:lag, :]
        elif lag > 0:
            a = norm_act[:-lag, :]
            b = norm_act[lag:, :]
        else:
            a = norm_act
            b = norm_act

        N = len(a)
        if N > 0:
            c = (a.T @ b) / max(1, N - 1)
            cross_corr[:, :, idx] = c

    return cross_corr, lags


def compute_shuffled_significance(activity, n_shuffles=100, min_shift_s=30.0, frame_dt=0.05, method='pearson', rng_seed=42):
    """
    Evalúa la significancia de la correlación de co-activación mediante desplazamientos temporales circulares.

    Parámetros
    ----------
    activity : ndarray, forma (n_frames, n_neurons)
    n_shuffles : int
    min_shift_s : float
        Shift mínimo en segundos.
    frame_dt : float
        Duración de frame en segundos.
    method : str
    rng_seed : int

    Devuelve
    --------
    results : dict
        - 'observed_corr': Matriz de correlación observada.
        - 'shuffled_mean': Media de la correlación nula.
        - 'shuffled_std': Desviación estándar de la correlación nula.
        - 'z_scores': (observada - media_nula) / std_nula.
        - 'p_values': p-valor empírico bilateral.
    """
    rng = np.random.default_rng(rng_seed)
    activity = np.asarray(activity, dtype=float)
    n_frames, n_neurons = activity.shape

    observed_corr = compute_pairwise_correlation(activity, method=method)

    min_shift = int(np.ceil(min_shift_s / frame_dt))
    max_shift = n_frames - min_shift

    if max_shift <= min_shift:
        # Sesion demasiado corta para shift de 30s
        min_shift = max(1, n_frames // 10)
        max_shift = n_frames - min_shift

    shuffled_corrs = np.zeros((n_shuffles, n_neurons, n_neurons), dtype=float)

    for s in range(n_shuffles):
        # Shift independiente por neurona
        shifted_activity = np.zeros_like(activity)
        for i in range(n_neurons):
            shift = rng.integers(min_shift, max_shift + 1) if max_shift > min_shift else rng.integers(1, n_frames)
            shifted_activity[:, i] = np.roll(activity[:, i], shift)

        shuffled_corrs[s] = compute_pairwise_correlation(shifted_activity, method=method)

    shuffled_mean = np.mean(shuffled_corrs, axis=0)
    shuffled_std = np.std(shuffled_corrs, axis=0)
    shuffled_std_safe = np.where(shuffled_std > 0, shuffled_std, 1.0)

    z_scores = (observed_corr - shuffled_mean) / shuffled_std_safe
    np.fill_diagonal(z_scores, 0.0)

    # p-valor empírico (cuántos shuffles igualan o superan la observada)
    greater_count = np.sum(shuffled_corrs >= observed_corr[None, :, :], axis=0)
    p_values = (greater_count + 1) / (n_shuffles + 1)
    np.fill_diagonal(p_values, 1.0)

    return {
        'observed_corr': observed_corr,
        'shuffled_mean': shuffled_mean,
        'shuffled_std': shuffled_std,
        'z_scores': z_scores,
        'p_values': p_values
    }


def find_coactive_assemblies(corr_matrix, n_clusters=None, distance_threshold=0.5):
    """
    Agrupa neuronas en ensambles co-activos usando clustering jerárquico.

    Parámetros
    ----------
    corr_matrix : ndarray, forma (n_neurons, n_neurons)
    n_clusters : int, opcional
    distance_threshold : float, opcional

    Devuelve
    --------
    order : ndarray
        Índices de neuronas reordenados por similitud.
    labels : ndarray
        Etiqueta de clúster para cada neurona.
    linkage_matrix : ndarray
        Matriz de enlace de scipy.cluster.hierarchy.
    """
    corr_clean = np.nan_to_num(corr_matrix, nan=0.0)
    # Convertir correlación en distancia d = 1 - r
    dist_matrix = np.clip(1.0 - corr_clean, 0.0, 2.0)
    np.fill_diagonal(dist_matrix, 0.0)

    # Matriz condensed para linkage
    condensed_dist = hierarchy.distance.squareform(dist_matrix, checks=False)
    linkage_matrix = hierarchy.linkage(condensed_dist, method='average')

    dendro = hierarchy.dendrogram(linkage_matrix, no_plot=True)
    order = np.array(dendro['leaves'], dtype=int)

    if n_clusters is not None:
        labels = hierarchy.fcluster(linkage_matrix, t=n_clusters, criterion='maxclust')
    else:
        labels = hierarchy.fcluster(linkage_matrix, t=distance_threshold, criterion='distance')

    return order, labels, linkage_matrix
