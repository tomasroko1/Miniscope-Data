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


def _find_mat_in_dir(folder):
    """Revisa si en una carpeta existe un archivo de trazas de CaliAli."""
    if not os.path.isdir(folder):
        return None
    try:
        entries = os.listdir(folder)
        # Prioridad 1: archivo que empiece con caliali_traces
        for fn in entries:
            fn_l = fn.lower()
            if fn_l.startswith("caliali_traces") and fn_l.endswith(".mat"):
                return os.path.join(folder, fn)
        # Prioridad 2: cualquier .mat que contenga 'traces' o 'trace'
        for fn in entries:
            fn_l = fn.lower()
            if fn_l.endswith(".mat") and ("traces" in fn_l or "trace" in fn_l):
                return os.path.join(folder, fn)
        # Prioridad 3: cualquier .mat de caliali que contenga 'result'
        for fn in entries:
            fn_l = fn.lower()
            if fn_l.endswith(".mat") and "caliali" in fn_l and "result" in fn_l:
                return os.path.join(folder, fn)
    except (OSError, IOError):
        pass
    return None


def find_caliali_traces_path(caliali_dir, animal, session_file, day_name=""):
    """Localiza el archivo CaliAli_traces_*.mat correspondiente a una sesion en GDi de forma instantanea."""
    if not caliali_dir or not os.path.exists(caliali_dir):
        return None

    date_match = re.search(r"\d{4}_\d{2}_\d{2}", session_file)
    date_str = date_match.group(0) if date_match else ""

    animal_clean = animal.split("_")[0].upper()
    try:
        animal_dirs = [
            d for d in os.listdir(caliali_dir)
            if os.path.isdir(os.path.join(caliali_dir, d))
            and (d.upper() == animal_clean or d.upper().startswith(animal_clean + "_") or animal_clean in d.upper())
        ]
    except (OSError, IOError):
        return None

    if not animal_dirs:
        return None

    clean_day = re.sub(r"[^A-Za-z0-9]", "", day_name).upper() if day_name else ""

    for ad in animal_dirs:
        ad_path = os.path.join(caliali_dir, ad)
        try:
            session_dirs = [d for d in os.listdir(ad_path) if os.path.isdir(os.path.join(ad_path, d))]
        except (OSError, IOError):
            continue

        # 1. Probar carpetas que coincidan con fecha o day_name
        candidates = []
        for sd in session_dirs:
            matches_date = bool(date_str and date_str in sd)
            clean_sd = re.sub(r"[^A-Za-z0-9]", "", sd).upper()
            matches_day = bool(clean_day and clean_day in clean_sd)
            if matches_date or matches_day:
                candidates.append(sd)

        # 2. Si no hubo coincidencia directa pero hay date_str, buscar en cualquier carpeta que tenga la fecha
        if not candidates and date_str:
            candidates = [sd for sd in session_dirs if date_str in sd]

        for sd in candidates:
            sd_path = os.path.join(ad_path, sd)
            sd_clean = sd.replace(" ", "_")

            # Patrón estándar universal según estructura de CaliAli:
            # CaliAli_traces_<ANIMAL_DIR>_<SESSION_DIR>.mat
            direct_patterns = [
                os.path.join(sd_path, "CaliAli_results", "CaliAli_traces_{}_{}.mat".format(ad, sd_clean)),
                os.path.join(sd_path, "CaliAli_results", "CaliAli_traces_{}_{}.mat".format(ad, sd)),
                os.path.join(sd_path, "caliali_results", "CaliAli_traces_{}_{}.mat".format(ad, sd_clean)),
                os.path.join(sd_path, "CaliAli_traces_{}_{}.mat".format(ad, sd_clean)),
                os.path.join(sd_path, "CaliAli_traces_{}_{}.mat".format(ad, sd)),
            ]
            for dp in direct_patterns:
                if os.path.isfile(dp):
                    return dp

            # Revisar carpetas típicas de resultados
            for subname in ("CaliAli_results", "caliali_results", "results", "output", "traces"):
                full_sub = os.path.join(sd_path, subname)
                found = _find_mat_in_dir(full_sub)
                if found:
                    return found

            # Revisar cualquier otra subcarpeta de la sesión (omitiendo bloques brutos de adquisición)
            try:
                for sub in os.listdir(sd_path):
                    if sub.upper() in ("A_AV1", "B_SAMPLE", "C_TEST", "D_AV2", "MINIAN", ".ZARR"):
                        continue
                    full_sub = os.path.join(sd_path, sub)
                    if os.path.isdir(full_sub):
                        found = _find_mat_in_dir(full_sub)
                        if found:
                            return found
            except (OSError, IOError):
                pass

            # Buscar directamente en la raíz de sd_path
            found = _find_mat_in_dir(sd_path)
            if found:
                return found

    return None


def _extract_sparse_or_dense(f, obj):
    """Extrae una matriz 2D (densa o sparse de MATLAB v7.3 HDF5) a numpy float32."""
    if isinstance(obj, h5py.Group):
        # Sparse matrix de MATLAB v7.3: tiene datasets 'data', 'ir', 'jc'
        if "data" in obj and "ir" in obj and "jc" in obj:
            data = np.array(obj["data"], dtype=np.float32)
            ir = np.array(obj["ir"], dtype=np.int64)
            jc = np.array(obj["jc"], dtype=np.int64)
            n_cols = len(jc) - 1
            n_rows = int(np.max(ir) + 1) if len(ir) > 0 else 0
            if "MATLAB_sparse" in obj.attrs:
                try:
                    n_rows = int(obj.attrs["MATLAB_sparse"])
                except Exception:
                    pass
            from scipy.sparse import csc_matrix
            mat = csc_matrix((data, ir, jc), shape=(n_rows, n_cols)).toarray()
            # En HDF5 de MATLAB, las matrices están transpuestas (filas x cols invertidas)
            return np.asarray(mat, dtype=np.float32)

        # Si el grupo tiene una clave interna 'S'
        if "S" in obj:
            return _extract_sparse_or_dense(f, obj["S"])

        raise ValueError("Grupo HDF5 no reconocido como matriz: {}".format(list(obj.keys())))

    elif isinstance(obj, h5py.Dataset):
        # Si es un dataset numérico regular
        if obj.dtype.kind in ("f", "i", "u"):
            return np.array(obj, dtype=np.float32)
        # Si es un dataset de referencias (Cell array en MATLAB)
        elif obj.dtype.kind == "O":
            mats = []
            for ref in np.asarray(obj).ravel():
                ref_obj = f[ref]
                mats.append(_extract_sparse_or_dense(f, ref_obj))
            return mats

    raise ValueError("Objeto no convertible a matriz: {}".format(type(obj)))


def _extract_session_vector(f, sess_ds):
    """Extrae un vector 1D de sesiones (numerico o strings) desde HDF5."""
    if isinstance(sess_ds, h5py.Dataset):
        if sess_ds.dtype.kind == "O":
            # Cell array de strings en MATLAB: cada celda es una referencia a un vector de caracteres
            labels = []
            for ref in np.asarray(sess_ds).ravel():
                char_ds = f[ref]
                s = "".join(chr(c) for c in np.array(char_ds).ravel())
                labels.append(s)
            return np.array(labels)
        elif sess_ds.dtype.kind in ("S", "U"):
            return np.array(sess_ds).ravel()
        else:
            return np.array(sess_ds).ravel()
    return np.array(sess_ds).ravel()


def read_caliali_traces(mat_path, signal_name="S"):
    """Lee CaliAli_traces_*.mat (formato MATLAB v7.3 HDF5 o v7 clásico).

    Args:
        mat_path: ruta al archivo .mat de CaliAli_traces
        signal_name: nombre del campo a extraer ('S', 'C_raw', 'C')

    Devuelve:
        subs_S: lista de arrays [sig_0, sig_1, sig_2, sig_3], cada uno (n_frames_i, n_neurons)
        fps: frecuencia de muestreo en Hz
    """
    fps = 20.0
    sig_key = signal_name.replace(".", "_")

    # Intento 1: Formato HDF5 (MATLAB v7.3)
    try:
        with h5py.File(mat_path, "r") as f:
            if "fps" in f:
                fps_ds = f["fps"]
                fps = float(fps_ds[0, 0] if fps_ds.ndim > 1 else fps_ds[0])
            elif "traces" in f and "fps" in f["traces"]:
                fps_ds = f["traces"]["fps"]
                fps = float(fps_ds[0, 0] if fps_ds.ndim > 1 else fps_ds[0])

            if "traces" not in f:
                raise KeyError("No se encontro el grupo 'traces' en {}. Claves: {}".format(mat_path, list(f.keys())))

            traces_obj = f["traces"]

            # Caso 1: traces es un Dataset de referencias (struct array de 4 subsesiones)
            if isinstance(traces_obj, h5py.Dataset) and traces_obj.dtype.kind == "O":
                subs_S = []
                for ref in np.asarray(traces_obj).ravel():
                    sub_grp = f[ref]
                    if sig_key not in sub_grp:
                        raise KeyError(f"Clave '{sig_key}' no encontrada en subsesion. Claves disponibles: {list(sub_grp.keys())}")
                    s_arr = _extract_sparse_or_dense(f, sub_grp[sig_key])
                    if s_arr.ndim == 2 and s_arr.shape[0] < s_arr.shape[1]:
                        s_arr = s_arr.T
                    subs_S.append(s_arr)
                return subs_S, fps

            # Caso 2: traces es un Grupo HDF5
            if isinstance(traces_obj, h5py.Group):
                # 2a. Sub-grupos explícitos (ej: 'sub1', 'sub2', ...)
                sub_keys = [k for k in traces_obj.keys() if k.isdigit() or k.lower().startswith("sub")]
                if sub_keys:
                    subs_S = []
                    for k in sorted(sub_keys):
                        s_arr = _extract_sparse_or_dense(f, traces_obj[k][sig_key])
                        if s_arr.ndim == 2 and s_arr.shape[0] < s_arr.shape[1]:
                            s_arr = s_arr.T
                        subs_S.append(s_arr)
                    return subs_S, fps

                # 2b. Matriz o Cell array con nombre sig_key
                if sig_key in traces_obj:
                    extracted = _extract_sparse_or_dense(f, traces_obj[sig_key])

                    # Si el campo era un cell array de 4 matrices, ya viene como lista
                    if isinstance(extracted, list):
                        subs_S = []
                        for sub_mat in extracted:
                            if sub_mat.ndim == 2 and sub_mat.shape[0] < sub_mat.shape[1]:
                                sub_mat = sub_mat.T
                            subs_S.append(sub_mat)
                        return subs_S, fps

                    s_arr = extracted
                    if s_arr.ndim == 2 and s_arr.shape[0] < s_arr.shape[1]:
                        s_arr = s_arr.T

                    # Particionar por traces.session si existe
                    if "session" in traces_obj:
                        sess_vec = _extract_session_vector(f, traces_obj["session"])
                        if len(sess_vec) == len(s_arr):
                            _, idx = np.unique(sess_vec, return_index=True)
                            unique_ordered = sess_vec[np.sort(idx)]
                            subs_S = [s_arr[sess_vec == u] for u in unique_ordered]
                            return subs_S, fps

                    return [s_arr], fps

            raise ValueError("Estructura no reconocida para 'traces' o campo '{}' en {}".format(sig_key, mat_path))

    except (OSError, IOError):
        # Intento 2: Formato MATLAB v5 / v7 clásico (scipy.io.loadmat)
        import scipy.io as sio
        d = sio.loadmat(mat_path, squeeze_me=True, struct_as_record=False)
        traces = d.get("traces")
        if traces is None:
            raise KeyError("No se encontro 'traces' en {}".format(mat_path))
        fps = float(getattr(d, "fps", getattr(traces, "fps", 20.0)))

        S = getattr(traces, sig_key)
        if hasattr(S, "toarray"):  # scipy sparse matrix
            S = S.toarray()

        if isinstance(S, np.ndarray) and S.dtype == object:
            subs_S = []
            for sub_s in S:
                arr = sub_s.toarray() if hasattr(sub_s, "toarray") else np.asarray(sub_s, dtype=np.float32)
                if arr.ndim == 2 and arr.shape[0] < arr.shape[1]:
                    arr = arr.T
                subs_S.append(arr)
            return subs_S, fps

        S = np.asarray(S, dtype=np.float32)
        if S.ndim == 2 and S.shape[0] < S.shape[1]:
            S = S.T

        sess = getattr(traces, "session", None)
        if sess is not None:
            sess_vec = np.asarray(sess).ravel()
            if len(sess_vec) == len(S):
                _, idx = np.unique(sess_vec, return_index=True)
                unique_ordered = sess_vec[np.sort(idx)]
                subs_S = [S[sess_vec == u] for u in unique_ordered]
                return subs_S, fps

        return [S], fps
