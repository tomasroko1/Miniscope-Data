import numpy as np
import matplotlib.pyplot as plt

def compute_rate_map(x, y, t, spikes, bin_size=2.5, sigma=None, arena_extent=None):
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
    
    # la arena circular (radio fijo de 25cm -> preguntar)
    if arena_extent is None:
        x_min, x_max, y_min, y_max = -25, 25, -25, 25
    else:
        x_min, x_max, y_min, y_max = arena_extent

    cx = (x_max + x_min) / 2
    cy = (y_max + y_min) / 2
    radio = (x_max - x_min) / 2

    xb = np.arange(x_min, x_max + bin_size, bin_size)
    yb = np.arange(y_min, y_max + bin_size, bin_size)
    xc = (xb[:-1] + xb[1:]) / 2
    yc = (yb[:-1] + yb[1:]) / 2
    XC, YC = np.meshgrid(xc, yc, indexing='ij')
    
    if sigma is None:
        sigma = bin_size

    # suavizado gaussiano antes de binear
    grid_x, grid_y = XC.ravel()[:, None], YC.ravel()[:, None]
    xd, yd = x[None, :], y[None, :]
    
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
    
    rate_map[dist_from_center > radio] = np.nan
    # tambien ocultamos pixeles donde el suavizado no llego (occ_s muy bajo)
    rate_map[occ_s < 1e-4] = np.nan

    return rate_map, xb, yb

def draw_trajectory(ax, x, y, spikes, neuron_id=""):
    ax.clear()
    idx = np.where(spikes > 0)[0]
    ax.plot(x, y, color='lightgray', linewidth=1, label='trayectoria')
    if len(idx) > 0:
        ax.scatter(x[idx], y[idx], color='red', s=12, zorder=5, label='spikes')
    
    ax.set_aspect('equal')
    ax.set_title(f"{neuron_id} - recorrido y spikes")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend(loc='upper right')

def draw_place_field(ax, rate_map, extent, neuron_id=""):
    ax.clear()
    m = np.nanmax(rate_map) if not np.isnan(rate_map).all() else 0.0
    cmap = plt.cm.jet.copy()
    cmap.set_bad(color='white')

    im = ax.imshow(rate_map.T, origin='lower', extent=extent, cmap=cmap, vmin=0, vmax=m if m > 0 else 1)
    ax.set_aspect('equal')
    ax.set_title(f"{neuron_id} - mapa (max: {m:.2f} hz)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    return im

def draw_neuron_analysis(ax1, ax2, x, y, spikes, rate_map, extent, neuron_id=""):
    draw_trajectory(ax1, x, y, spikes, neuron_id)
    return draw_place_field(ax2, rate_map, extent, neuron_id)
