import os
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor

from core.data_loader import list_sessions, load_session
from core.place_fields import compute_rate_map, draw_neuron_analysis

# =============================================================================
# PARAMETROS
# =============================================================================
BIN_SIZE = 2.5
SIGMA = 2.0
OUTPUT_DIR = "results"
MAX_WORKERS = 4 
SAVE_PLOTS = True
# =============================================================================

def procesar_sesion(session_info):
    animal = session_info['animal']
    session_file = session_info['session_file']
    
    print(f"[{animal}] Procesando {session_file}...")
    
    # cargar datos
    data = load_session(animal, session_file)
    
    # crear path
    cache_dir = os.path.join(OUTPUT_DIR, "cache", animal)
    os.makedirs(cache_dir, exist_ok=True)
    
    if SAVE_PLOTS:
        plot_dir = os.path.join(OUTPUT_DIR, "plots", animal, session_file.replace(".mat", ""))
        os.makedirs(plot_dir, exist_ok=True)
    
    resultados_sesion = {
        'animal': animal,
        'session_file': session_file,
        'bin_size': BIN_SIZE,
        'sigma': SIGMA,
        'rate_maps': {}
    }
    
    # manejar sesiones simples vs merged
    subsessions = data['subsessions'] if data['type'] == 'merged' else [data]
    
    for sub_idx, sub in enumerate(subsessions):
        x, y, t, spikes = sub['x'], sub['y'], sub['t'], sub['S']
        n_neurons = sub['n_neurons']
        
        # guardamos por neurona
        for n in range(n_neurons):
            neuron_spikes = spikes[:, n]
            if np.sum(neuron_spikes) == 0:
                continue
                
            rate_map, x_bins, y_bins = compute_rate_map(
                x, y, t, neuron_spikes,
                bin_size=BIN_SIZE,
                sigma=SIGMA
            )
            extent = [x_bins[0], x_bins[-1], y_bins[0], y_bins[-1]]
            
            # guardamos el mapa en la memoria de resultados
            key = f"sub_{sub_idx}_n_{n}" if data['type'] == 'merged' else f"n_{n}"
            resultados_sesion['rate_maps'][key] = rate_map
            
            if SAVE_PLOTS:
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
                title = f"Neurona {n}" + (f" (Sub {sub_idx})" if data['type'] == 'merged' else "")
                draw_neuron_analysis(ax1, ax2, x, y, neuron_spikes, rate_map, extent, title)
                fig.savefig(os.path.join(plot_dir, f"{key}.png"), dpi=100, bbox_inches='tight')
                plt.close(fig)
    
    # exportar el .pkl unificado de la sesion
    out_pkl = os.path.join(cache_dir, session_file.replace(".mat", ".pkl"))
    with open(out_pkl, 'wb') as f:
        pickle.dump(resultados_sesion, f)
        
    print(f"[{animal}] Terminando {session_file}. Guardado en {out_pkl}")

def main():
    print(f"Iniciando procesamiento...")
    sessions = list_sessions()
    
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        executor.map(procesar_sesion, sessions)
        
    print("Batch finalizado.")

if __name__ == "__main__":
    main()
