"""caliali_loader.py

Modulo para localizar y cargar datos procesados por CaliAli en GDi/
(archivos MATLAB v7.3 HDF5 CaliAli_traces_*.mat) y alinearlos con el tracking.
"""

from __future__ import print_function
import os
import re
from pathlib import Path
import h5py
import numpy as np


def find_caliali_traces_path(caliali_dir, animal, session_file, day_name=""):
    """Localiza el archivo CaliAli_traces_*.mat correspondiente a una sesion en GDi."""
    if not caliali_dir or not os.path.exists(caliali_dir):
        return None

    date_match = re.search(r"\d{4}_\d{2}_\d{2}", session_file)
    date_str = date_match.group(0) if date_match else ""

    animal_clean = animal.split("_")[0].upper()
    animal_dirs = [
        d for d in os.listdir(caliali_dir)
        if os.path.isdir(os.path.join(caliali_dir, d))
        and (d.upper() == animal_clean or d.upper().startswith(animal_clean + "_"))
    ]

    if not animal_dirs:
        return None

    for ad in animal_dirs:
        ad_path = os.path.join(caliali_dir, ad)
        for sd in os.listdir(ad_path):
            sd_path = os.path.join(ad_path, sd)
            if not os.path.isdir(sd_path):
                continue
            matches_date = bool(date_str and date_str in sd)
            matches_day = bool(day_name and day_name.replace("_", "").upper() in sd.replace("_", "").upper())
            if matches_date or matches_day:
                results_dir = os.path.join(sd_path, "CaliAli_results")
                if os.path.isdir(results_dir):
                    for fn in os.listdir(results_dir):
                        if fn.startswith("CaliAli_traces_") and fn.endswith(".mat"):
                            return os.path.join(results_dir, fn)
    return None


def read_caliali_traces(mat_path):
    """Lee CaliAli_traces_*.mat (formato MATLAB v7.3 HDF5).

    Devuelve:
        subs_S: lista de arrays [S_0, S_1, S_2, S_3], cada uno (n_frames_i, n_neurons)
        fps: frecuencia de muestreo en Hz
    """
    with h5py.File(mat_path, "r") as f:
        fps = 20.0
        if "fps" in f:
            fps_ds = f["fps"]
            fps = float(fps_ds[0, 0] if fps_ds.ndim > 1 else fps_ds[0])
        elif "traces" in f and "fps" in f["traces"]:
            fps_ds = f["traces"]["fps"]
            fps = float(fps_ds[0, 0] if fps_ds.ndim > 1 else fps_ds[0])

        if "traces" not in f:
            raise KeyError("No se encontro el grupo 'traces' en {}. Claves: {}".format(mat_path, list(f.keys())))

        traces_obj = f["traces"]

        # Caso 1: struct array de referencias (4 subsesiones)
        if isinstance(traces_obj, h5py.Dataset) and traces_obj.dtype.kind == "O":
            subs_S = []
            for ref in traces_obj.flat:
                sub_grp = f[ref]
                s_ds = sub_grp["S"]
                s_arr = np.array(s_ds, dtype=np.float32)
                if s_arr.ndim == 2 and s_arr.shape[0] < s_arr.shape[1]:
                    s_arr = s_arr.T
                subs_S.append(s_arr)
            return subs_S, fps

        # Caso 2: Grupo traces
        if isinstance(traces_obj, h5py.Group):
            sub_keys = [k for k in traces_obj.keys() if k.isdigit() or k.lower().startswith("sub")]
            if sub_keys:
                subs_S = []
                for k in sorted(sub_keys):
                    s_arr = np.array(traces_obj[k]["S"], dtype=np.float32)
                    if s_arr.ndim == 2 and s_arr.shape[0] < s_arr.shape[1]:
                        s_arr = s_arr.T
                    subs_S.append(s_arr)
                return subs_S, fps

            # Matriz 2D global 'S' particionada por 'session'
            if "S" in traces_obj:
                s_arr = np.array(traces_obj["S"], dtype=np.float32)
                if s_arr.ndim == 2 and s_arr.shape[0] < s_arr.shape[1]:
                    s_arr = s_arr.T

                if "session" in traces_obj:
                    sess_ds = traces_obj["session"]
                    sess_vec = np.array(sess_ds).ravel()
                    if len(sess_vec) == len(s_arr):
                        _, idx = np.unique(sess_vec, return_index=True)
                        unique_ordered = sess_vec[np.sort(idx)]
                        subs_S = [s_arr[sess_vec == u] for u in unique_ordered]
                        return subs_S, fps

                return [s_arr], fps

        raise ValueError("Estructura no reconocida para 'traces' en {}".format(mat_path))
