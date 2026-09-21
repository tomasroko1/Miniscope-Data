import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.ndimage import gaussian_filter


def threshold_event_frames(signal, threshold_source=None, threshold=None,
                           threshold_sd=3.0):
    """Convierte S deconvolucionada en frames-evento binarios.

    Un frame vale 1 cuando ``S > threshold_sd * SD(S)`` y 0 en caso
    contrario. La SD se calcula por celula y sesion.
    """
    signal = np.asarray(signal, dtype=float)
    clean = np.where(np.isfinite(signal), signal, 0.0)
    if threshold is None:
        source = signal if threshold_source is None else np.asarray(
            threshold_source, dtype=float
        )
        source = source[np.isfinite(source)]
        if source.size < 2:
            raise ValueError("insuficientes datos")
        threshold = threshold_sd * float(np.std(source, ddof=1))
    threshold = float(threshold)
    events = (clean > threshold).astype(float)
    return events, threshold


def event_run_lengths(events):
    """Longitud en frames de cada racha contigua de frames-evento."""
    mask = np.asarray(events, dtype=float) > 0
    starts = np.flatnonzero(mask & ~np.r_[False, mask[:-1]])
    ends = np.flatnonzero(mask & ~np.r_[mask[1:], False])
    return ends - starts + 1

def _validate_event_frames(events, expected_length):
    events = np.asarray(events, dtype=float)
    if events.ndim != 1 or len(events) != expected_length:
        raise ValueError("events debe tener un valor por frame de posicion")
    finite = events[np.isfinite(events)]
    if not np.all((finite == 0) | (finite == 1)):
        raise ValueError(
            "draw_trajectory espera frames-evento binarios; "
            "use threshold_event_frames(S) antes de graficar"
        )
    return events


def draw_trajectory(ax, x, y, events, neuron_id=""):
    """Dibuja trayectoria y frames-evento."""
    ax.clear()
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y):
        raise ValueError("x e y deben ser vectores de igual longitud")
    events = _validate_event_frames(events, len(x))
    valid_position = np.isfinite(x) & np.isfinite(y)
    event_position = valid_position & (events > 0)

    ax.plot(x[valid_position], y[valid_position], color="lightgray",
            linewidth=1, label="trayectoria")
    if np.any(event_position):
        ax.scatter(x[event_position], y[event_position], color="red", s=12,
                   zorder=5, label="frames-evento")
    ax.set_aspect("equal")
    ax.set_title(f"{neuron_id} - trayectoria y eventos")
    ax.set_xlabel("x [cm]")
    ax.set_ylabel("y [cm]")
    ax.legend(loc="upper right")


def compute_rate_map(x, y, t, spikes, bin_size=1, sigma=1.5):
    # calcula la duración de cada frame
    dt = np.zeros(len(t))
    if len(t) > 1:
        dt[:-1] = np.diff(t)
        dt[-1] = dt[-2]
    else:
        dt[:] = 0.05

    # frames validos. dropeo si hay nulos en x,y,t o spikes
    valid = ~(np.isnan(x) | np.isnan(y) | np.isnan(t) | np.isnan(spikes))
    x, y, dt, spikes = x[valid], y[valid], dt[valid], spikes[valid]

    # la arena es circular (radio fijo de 25cm)
    x_min, x_max, y_min, y_max = -25, 25, -25, 25

    cx = (x_max + x_min) / 2 # cx es el centro de la arena
    cy = (y_max + y_min) / 2 # cy es el centro de la arena
    radio = (x_max - x_min) / 2

    xb = np.arange(x_min, x_max + bin_size, bin_size) # armo un grid de bines de x. x_min es el borde izquierdo del primer bin, x_max es el borde derecho del ultimo bin, le sumo bin_size para que el ultimo bin tenga el borde derecho en x_max
    yb = np.arange(y_min, y_max + bin_size, bin_size) # bines de y
    xc = (xb[:-1] + xb[1:]) / 2 # xc es el centro de cada bin de x
    yc = (yb[:-1] + yb[1:]) / 2 # yc es el centro de cada bin de y
    XC, YC = np.meshgrid(xc, yc, indexing='ij') # XC y YC son matrices que contienen las coordenadas de cada bin en el grid. XC tiene las coordenadas x de cada bin, YC tiene las coordenadas y de cada bin. el indexing='ij' es para que el primer indice sea el de x y el segundo el de y.

    ## suavizado gaussiano antes de binear

    grid_x, grid_y = XC.ravel()[:, None], YC.ravel()[:, None] # grid_x y grid_y son vectores columna con las coordenadas de cada bin en el grid. el ravel() es para aplanar la matriz, y el [:, None] es para convertirlo en un vector columna.
    xd, yd = x[None, :], y[None, :] # xd y yd son vectores fila con las coordenadas de cada frame. el [None, :] es para convertirlo en un vector fila.

    # proyectamos cada dato como una campana de gauss al centro de cada cuadradito
    gauss = np.exp(-((grid_x - xd)**2 + (grid_y - yd)**2) / (2 * sigma**2))
    occ_s = (gauss @ dt).reshape(XC.shape)
    spk_s = (gauss @ spikes).reshape(XC.shape)

    rate_map = np.zeros_like(occ_s, dtype=float)
    mask = occ_s > 1e-4
    rate_map[mask] = spk_s[mask] / occ_s[mask]

    # armamos el circulo de bines para el plot
    xc = (xb[:-1] + xb[1:]) / 2
    yc = (yb[:-1] + yb[1:]) / 2
    XC, YC = np.meshgrid(xc, yc, indexing='ij')
    dist_from_center = np.sqrt((XC - cx)**2 + (YC - cy)**2)
    
    # oculto los pixeles de la matriz cuadrada que quedan fuera de la arena circular
    rate_map[dist_from_center > radio] = np.nan
    # tambien ocultamos pixeles donde el suavizado no llego (occ_s muy bajo)
    rate_map[occ_s < 1e-4] = np.nan

    return rate_map, xb, yb


def gaussian_rate_maps_matlab(x, y, t, activity, bin_size=1.0,
                              sigma=5.0, arena_radius=25.0,
                              block_size=256):
    """Rate maps from the continuous-position Gaussian estimator.

    This is the direct vectorised equivalent of the MATLAB calculation::

        W = exp(-(xdiff + ydiff) / sigma_sq);
        coverage = sum(W, 2);
        mapr = sample_rate * (W * r) ./ coverage;

    Here ``sigma`` is the usual Gaussian standard deviation in cm, hence
    ``sigma_sq`` in the expression above is ``2 * sigma**2``.  The sampling
    rate is obtained from the timestamps rather than assumed: the HabL data
    are sampled at 20 Hz, not 50 Hz.  This only changes the units of the map,
    not its spatial pattern.

    ``activity`` is one value per frame and cell, normally a binary event
    matrix.  Evaluation pixels are processed in blocks, so the result is
    mathematically the same as allocating the very large full ``W`` matrix
    but does not require hundreds of MB of temporary memory.

    Returns
    -------
    maps : ndarray, shape (n_cells, n_x_bins, n_y_bins)
        Activity rate in events/s (or activity units/s).
    edges : ndarray
        Shared x/y grid edges in cm.
    metadata : dict
        Includes the sampling rate, the literal MATLAB ``sigma_sq`` and
        coverage map for provenance.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    t = np.asarray(t, dtype=float)
    activity = np.asarray(activity, dtype=float)
    if activity.ndim == 1:
        activity = activity[:, None]
    if activity.ndim != 2 or activity.shape[0] != len(x):
        raise ValueError("activity debe tener una fila por frame")
    if sigma <= 0 or bin_size <= 0 or block_size <= 0:
        raise ValueError("sigma, bin_size y block_size deben ser positivos")

    # Frames absent from an imaging block are absent for every local cell in
    # these data.  Dropping them is preferable to treating them as no-event.
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(t)
    valid &= np.isfinite(activity).all(axis=1)
    x, y, t, activity = x[valid], y[valid], t[valid], activity[valid]
    if len(t) < 2:
        raise ValueError("no hay suficientes frames validos")
    frame_dt = np.diff(t)
    positive_dt = frame_dt[frame_dt > 0]
    if positive_dt.size == 0:
        raise ValueError("los timestamps no son estrictamente crecientes")
    sample_rate_hz = 1.0 / float(np.median(positive_dt))

    edges = np.arange(-arena_radius, arena_radius + bin_size, bin_size)
    centres = (edges[:-1] + edges[1:]) / 2.0
    xx, yy = np.meshgrid(centres, centres, indexing="ij")
    pixels = np.column_stack((xx.ravel(), yy.ravel())).astype(np.float32)
    positions = np.column_stack((x, y)).astype(np.float32)
    values = activity.astype(np.float32, copy=False)
    sigma_sq = np.float32(2.0 * sigma ** 2)

    n_pixels = len(pixels)
    n_cells = values.shape[1]
    numerators = np.empty((n_pixels, n_cells), dtype=np.float32)
    coverage = np.empty(n_pixels, dtype=np.float32)
    for start in range(0, n_pixels, block_size):
        stop = min(start + block_size, n_pixels)
        dx = pixels[start:stop, 0, None] - positions[None, :, 0]
        dy = pixels[start:stop, 1, None] - positions[None, :, 1]
        weights = np.exp(-(dx * dx + dy * dy) / sigma_sq).astype(
            np.float32, copy=False
        )
        coverage[start:stop] = weights.sum(axis=1)
        numerators[start:stop] = weights @ values

    rates = np.divide(
        sample_rate_hz * numerators,
        coverage[:, None],
        out=np.full_like(numerators, np.nan),
        where=coverage[:, None] > 1e-4,
    ).T.reshape(n_cells, len(centres), len(centres))
    arena = xx ** 2 + yy ** 2 <= arena_radius ** 2
    rates[:, ~arena] = np.nan
    return rates, edges, {
        "sample_rate_hz": sample_rate_hz,
        "sigma_cm": float(sigma),
        "sigma_sq_cm2": float(sigma_sq),
        "coverage": coverage.reshape(xx.shape),
        "valid_frames": int(valid.sum()),
        "excluded_frames": int((~valid).sum()),
    }


def draw_place_field(ax, rate_map, extent, neuron_id=""):
    """Dibuja un event-rate map ya calculado, en eventos por segundo."""
    ax.clear()
    m = np.nanmax(rate_map) if not np.isnan(rate_map).all() else 0.0
    cmap = plt.cm.turbo.copy()
    cmap.set_bad(color="white")
    im = ax.imshow(
        rate_map.T, origin="lower", extent=extent, cmap=cmap,
        vmin=0, vmax=m if m > 0 else 1,
    )
    ax.set_aspect("equal")
    ax.set_title(f"{neuron_id} - event-rate map (max: {m:.2f} eventos/s)")
    ax.set_xlabel("x [cm]")
    ax.set_ylabel("y [cm]")
    return im


def draw_neuron_analysis(ax1, ax2, x, y, events, rate_map, extent,
                         neuron_id=""):
    """Dibuja trayectoria/eventos y su event-rate map asociado."""
    draw_trajectory(ax1, x, y, events, neuron_id)
    return draw_place_field(ax2, rate_map, extent, neuron_id)


def make_rate_grid(grid_size=1.0, arena_radius=25.0, margin=10.0):
    edges = np.arange(-arena_radius, arena_radius + grid_size, grid_size)
    centres = (edges[:-1] + edges[1:]) / 2
    xx, yy = np.meshgrid(centres, centres, indexing="ij")
    arena = xx ** 2 + yy ** 2 <= arena_radius ** 2

    full_edges = np.arange(
        -arena_radius - margin,
        arena_radius + margin + grid_size,
        grid_size,
    )
    full_centres = (full_edges[:-1] + full_edges[1:]) / 2
    crop = np.flatnonzero(
        (full_centres >= -arena_radius) & (full_centres < arena_radius)
    )
    return dict(edges=edges, full_edges=full_edges, crop=crop, arena=arena,
                grid_size=grid_size)


def _crop_map(values, crop):
    return values[np.ix_(crop, crop)]


def _smooth_occupancy(occupancy, prepared):
    smoothed = gaussian_filter(occupancy, prepared["sigma_bins"], mode="constant")
    smoothed = _crop_map(smoothed, prepared["crop"])
    smoothed[~prepared["arena"]] = 0.0
    return smoothed


def prepare_rate_session(sub, grid_size=1.0, sigma=1.5,
                         arena_radius=25.0, margin=10.0):
    grid = make_rate_grid(grid_size, arena_radius, margin)
    full_edges = grid["full_edges"]
    x = np.asarray(sub["x"], dtype=float)
    y = np.asarray(sub["y"], dtype=float)
    t = np.asarray(sub["t"], dtype=float)
    S = np.asarray(sub["S"], dtype=float)

    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(t)
    valid &= np.isfinite(S).all(axis=1)
    x, y, t, S = x[valid], y[valid], t[valid], S[valid]
    if len(t) < 2:
        raise ValueError("no hay suficientes frames validos")

    dt = np.empty(len(t), dtype=float)
    dt[:-1] = np.diff(t)
    positive_dt = dt[:-1][dt[:-1] > 0]
    dt[-1] = np.median(positive_dt) if positive_dt.size else 0.05
    if np.any(dt <= 0):
        raise ValueError("los timestamps no son estrictamente crecientes")

    n_side = len(full_edges) - 1
    first_centre = full_edges[0] + grid_size / 2
    ux = (x - first_centre) / grid_size
    uy = (y - first_centre) / grid_size
    ix = np.floor(ux).astype(int)
    iy = np.floor(uy).astype(int)
    inside = (ix >= 0) & (ix + 1 < n_side) & (iy >= 0) & (iy + 1 < n_side)
    x, y, t, dt, S, ix, iy, ux, uy = (
        value[inside] for value in (x, y, t, dt, S, ix, iy, ux, uy)
    )

    wx = ux - ix
    wy = uy - iy
    rows = np.concatenate((
        ix * n_side + iy,
        (ix + 1) * n_side + iy,
        ix * n_side + iy + 1,
        (ix + 1) * n_side + iy + 1,
    ))
    columns = np.tile(np.arange(len(t)), 4)
    weights = np.concatenate((
        (1 - wx) * (1 - wy),
        wx * (1 - wy),
        (1 - wx) * wy,
        wx * wy,
    ))
    frame_grid = sparse.csr_matrix(
        (weights, (rows, columns)),
        shape=(n_side * n_side, len(t)),
    )

    midpoint = t[0] + (t[-1] - t[0]) / 2
    first = t <= midpoint
    second = ~first
    occupancy = np.asarray(frame_grid @ dt).reshape(n_side, n_side)
    occ_first = np.asarray(frame_grid[:, first] @ dt[first]).reshape(n_side, n_side)
    occ_second = np.asarray(frame_grid[:, second] @ dt[second]).reshape(n_side, n_side)

    prepared = dict(
        **grid, S=S, x=x, y=y, t=t, dt=dt, frame_grid=frame_grid,
        first=first, second=second, full_shape=occupancy.shape,
        sigma_bins=sigma / grid_size,
    )
    prepared["occupancy"] = _smooth_occupancy(occupancy, prepared)
    prepared["occ_first"] = _smooth_occupancy(occ_first, prepared)
    prepared["occ_second"] = _smooth_occupancy(occ_second, prepared)
    return prepared


def rate_map_from_activity(activity_flat, occupancy, prepared):
    activity = np.asarray(activity_flat).reshape(prepared["full_shape"])
    numerator = gaussian_filter(activity, prepared["sigma_bins"], mode="constant")
    numerator = _crop_map(numerator, prepared["crop"])
    rate = np.full(prepared["arena"].shape, np.nan)
    ok = prepared["arena"] & (occupancy > 1e-4)
    rate[ok] = numerator[ok] / occupancy[ok]
    return rate


def rate_maps_from_activity(activity_maps, occupancy, prepared):
    numerator = gaussian_filter(
        activity_maps,
        (0, prepared["sigma_bins"], prepared["sigma_bins"]),
        mode="constant",
    )
    crop = prepared["crop"]
    numerator = numerator[:, crop][:, :, crop]
    rates = np.full_like(numerator, np.nan)
    ok = prepared["arena"] & (occupancy > 1e-4)
    rates[:, ok] = numerator[:, ok] / occupancy[ok]
    return rates


def spatial_information(rate, occupancy, arena):
    ok = arena & (occupancy > 0) & np.isfinite(rate)
    if ok.sum() < 2 or occupancy[ok].sum() <= 0:
        return np.nan
    p = occupancy[ok] / occupancy[ok].sum()
    r = rate[ok]
    mean_rate = np.sum(p * r)
    if mean_rate <= 0:
        return 0.0
    ratio = r / mean_rate
    positive = ratio > 0
    return float(np.sum(p[positive] * ratio[positive] * np.log2(ratio[positive])))


def spatial_information_stack(rates, occupancy, arena):
    ok = arena & (occupancy > 1e-4)
    p = occupancy[ok] / occupancy[ok].sum()
    values = rates[:, ok]
    means = np.sum(values * p[None, :], axis=1)
    ratio = np.divide(values, means[:, None], out=np.zeros_like(values),
                      where=means[:, None] > 0)
    terms = np.zeros_like(ratio)
    positive = ratio > 0
    terms[positive] = ratio[positive] * np.log2(ratio[positive])
    return np.sum(terms * p[None, :], axis=1)


def map_correlation(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return np.nan
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def spatial_maps_and_metrics(signal, prepared):
    """Mapas completo/mitades y metricas para una actividad ya definida.

    En el pipeline de imaging ``signal`` debe ser la salida binaria de
    :func:`threshold_event_frames`.
    """
    signal = np.asarray(signal, dtype=float)
    frame_grid = prepared["frame_grid"]
    first = prepared["first"]
    second = prepared["second"]
    whole = rate_map_from_activity(
        frame_grid @ signal, prepared["occupancy"], prepared
    )
    first_map = rate_map_from_activity(
        frame_grid[:, first] @ signal[first], prepared["occ_first"], prepared
    )
    second_map = rate_map_from_activity(
        frame_grid[:, second] @ signal[second], prepared["occ_second"], prepared
    )
    return {
        "whole_map": whole,
        "first_map": first_map,
        "second_map": second_map,
        "information": spatial_information(
            whole, prepared["occupancy"], prepared["arena"]
        ),
        "stability": map_correlation(first_map, second_map),
        "peak": float(np.nanmax(whole)),
        "activity_sum": float(np.sum(signal)),
        "nonzero_frames": int(np.count_nonzero(signal)),
    }


def _rowwise_correlation(a, b):
    ac = a - a.mean(axis=1, keepdims=True)
    bc = b - b.mean(axis=1, keepdims=True)
    denominator = np.sqrt(np.sum(ac * ac, axis=1) * np.sum(bc * bc, axis=1))
    return np.divide(np.sum(ac * bc, axis=1), denominator,
                     out=np.full(len(a), np.nan), where=denominator > 0)


def analyse_spatial_cell(signal, prepared, rng, n_shuffles=100, min_shift_s=30.0):
    frame_grid = prepared["frame_grid"]
    dt = prepared["dt"]
    first = prepared["first"]
    second = prepared["second"]
    occupancy = prepared["occupancy"]
    occ_first = prepared["occ_first"]
    occ_second = prepared["occ_second"]

    observed = spatial_maps_and_metrics(signal, prepared)
    whole_map = observed["whole_map"]
    first_map = observed["first_map"]
    second_map = observed["second_map"]
    information = observed["information"]
    stability = observed["stability"]

    frame_dt = np.median(dt)
    min_shift = int(np.ceil(min_shift_s / frame_dt))
    max_shift = len(signal) - min_shift
    if max_shift <= min_shift:
        raise ValueError("sesion muy corta")
    shifts = rng.integers(min_shift, max_shift + 1, size=n_shuffles)
    indices = (np.arange(len(signal))[:, None] - shifts[None, :]) % len(signal)
    shifted = signal[indices]

    shape = prepared["full_shape"]
    whole_activity = np.asarray(frame_grid @ shifted).T.reshape(n_shuffles, *shape)
    first_activity = np.asarray(frame_grid[:, first] @ shifted[first]).T.reshape(
        n_shuffles, *shape
    )
    second_activity = np.asarray(frame_grid[:, second] @ shifted[second]).T.reshape(
        n_shuffles, *shape
    )
    shuffled_whole = rate_maps_from_activity(whole_activity, occupancy, prepared)
    shuffled_first = rate_maps_from_activity(first_activity, occ_first, prepared)
    shuffled_second = rate_maps_from_activity(second_activity, occ_second, prepared)
    shuffled_information = spatial_information_stack(
        shuffled_whole, occupancy, prepared["arena"]
    )
    common = prepared["arena"] & (occ_first > 1e-4) & (occ_second > 1e-4)
    shuffled_stability = _rowwise_correlation(
        shuffled_first[:, common], shuffled_second[:, common]
    )
    return dict(
        whole_map=whole_map,
        first_map=first_map,
        second_map=second_map,
        information=information,
        stability=stability,
        shuffled_information=shuffled_information,
        shuffled_stability=shuffled_stability,
        activity_sum=float(np.sum(signal)),
    )


def shuffled_stability_distribution(signal, prepared, rng, n_shuffles=100,
                                     min_shift_s=30.0):
    """Nulo de estabilidad para una senal sin recalcular informacion espacial."""
    frame_grid = prepared["frame_grid"]
    dt = prepared["dt"]
    first = prepared["first"]
    second = prepared["second"]

    frame_dt = np.median(dt)
    min_shift = int(np.ceil(min_shift_s / frame_dt))
    max_shift = len(signal) - min_shift
    if max_shift <= min_shift:
        raise ValueError("sesion muy corta")
    shifts = rng.integers(min_shift, max_shift + 1, size=n_shuffles)
    indices = (np.arange(len(signal))[:, None] - shifts[None, :]) % len(signal)
    shifted = signal[indices]

    shape = prepared["full_shape"]
    first_activity = np.asarray(
        frame_grid[:, first] @ shifted[first]
    ).T.reshape(n_shuffles, *shape)
    second_activity = np.asarray(
        frame_grid[:, second] @ shifted[second]
    ).T.reshape(n_shuffles, *shape)
    shuffled_first = rate_maps_from_activity(
        first_activity, prepared["occ_first"], prepared
    )
    shuffled_second = rate_maps_from_activity(
        second_activity, prepared["occ_second"], prepared
    )
    common = (
        prepared["arena"]
        & (prepared["occ_first"] > 1e-4)
        & (prepared["occ_second"] > 1e-4)
    )
    return _rowwise_correlation(
        shuffled_first[:, common], shuffled_second[:, common]
    )
