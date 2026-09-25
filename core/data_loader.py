import os
import scipy.io
import numpy as np

def get_data_dir():
    configured_dir = os.environ.get('MINISCOPE_DATA_DIR')
    if configured_dir:
        return os.path.abspath(os.path.expanduser(configured_dir))

    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_dir, 'data')

def list_sessions(data_dir=None):
    # lista todas las sesiones disponibles devolviendo una lista de diccionarios
    if data_dir is None:
        data_dir = get_data_dir()
    
    sessions = []
    if not os.path.exists(data_dir):
        return sessions
        
    for animal in sorted(os.listdir(data_dir)):
        animal_path = os.path.join(data_dir, animal)
        if not os.path.isdir(animal_path):
            continue
            
        for file_name in sorted(os.listdir(animal_path)):
            if not file_name.endswith(".mat"):
                continue
            sessions.append({
                'animal': animal,
                'session_file': file_name
            })
    return sessions

def load_session(animal, session_file, data_dir=None):
    # carga el .mat file y devuelve un diccionario con los datos ordenados de la session.
    if data_dir is None:
        data_dir = get_data_dir()

    fpath = os.path.join(data_dir, animal, session_file)
    if not os.path.exists(fpath):
        raise FileNotFoundError(f"Session not found at {fpath}")

    d = scipy.io.loadmat(fpath, squeeze_me=True, struct_as_record=False)
    if 'pos' not in d or 'act' not in d:
        raise ValueError(f"File {session_file} does not contain pos/act fields")

    pos = d['pos']
    act = d['act']
    sess = d.get('sess', None)

    x_raw = np.asarray(pos.x)
    is_merged = x_raw.dtype == object or (x_raw.ndim == 1 and isinstance(x_raw[0], np.ndarray))

    if not is_merged:
        # simple session
        x = np.asarray(pos.x, dtype=float)
        y = np.asarray(pos.y, dtype=float)
        v = np.asarray(pos.v, dtype=float) if hasattr(pos, 'v') else None
        
        ## TODO chequear este 0.05, me depende de los frames del registro!
        t = np.asarray(act.t, dtype=float) if hasattr(act, 't') else np.arange(len(x)) * 0.05
        S = np.asarray(act.S, dtype=float)
        
        if S.ndim == 1:
            S = S.reshape(-1, 1)

        return {
            'type': 'simple',
            'animal': animal,
            'session_file': session_file,
            'x': x,
            'y': y,
            'v': v,
            't': t,
            'S': S,
            'n_neurons': S.shape[1],
            'sess': sess
        }
    else:
        # multi-subsession merged file
        n_subs = len(x_raw)
        mapping = np.asarray(act.mapping) if hasattr(act, 'mapping') else None
        
        subs_data = []
        for i in range(n_subs):
            sx = np.asarray(pos.x[i], dtype=float)
            sy = np.asarray(pos.y[i], dtype=float)
            sv = np.asarray(pos.v[i], dtype=float) if hasattr(pos, 'v') else None
            st = np.asarray(act.t[i], dtype=float) if hasattr(act, 't') else np.arange(len(sx)) * 0.05
            sS = np.asarray(act.S[i], dtype=float)
            if sS.ndim == 1:
                sS = sS.reshape(-1, 1)

            subs_data.append({
                'x': sx, 'y': sy, 'v': sv, 't': st, 'S': sS,
                'n_neurons': sS.shape[1]
            })

        return {
            'type': 'merged',
            'animal': animal,
            'session_file': session_file,
            'subsessions': subs_data,
            'mapping': mapping,
            'n_subsessions': n_subs,
            'sess': sess
        }


def get_common_neurons(mapping):
    # devuelve los indices de neuronas que estan presente en todas las sub sessiones.
    if mapping is None:
        return np.array([], dtype=int)
    mapping = np.asarray(mapping)
    valid_mask = ~np.isnan(mapping).any(axis=1)
    return np.where(valid_mask)[0]


def mapping_rows_for_columns(mapping, sub_index, n_columns):
    """devuelve la fila global correspondiente a cada columna de s."""
    mapping = np.asarray(mapping, dtype=float)
    local_ids = mapping[:, sub_index]
    rows = np.flatnonzero(np.isfinite(local_ids))
    if len(rows) != n_columns:
        raise ValueError(
            f"mapping ambiguo en la subsesion {sub_index + 1}: "
            f"hay {len(rows)} ids locales y {n_columns} columnas en s"
        )
    ordered_rows = rows[np.argsort(local_ids[rows])]
    ordered_ids = local_ids[ordered_rows]
    if len(np.unique(ordered_ids)) != len(ordered_ids):
        raise ValueError(f"hay ids locales duplicados en la subsesion {sub_index + 1}")
    return ordered_rows.astype(int)


def mapping_columns(mapping, global_ids, sub_index, n_columns):
    """Translate local mapping IDs to columns of a filtered activity matrix.

    In this dataset ``mapping`` stores zero-based IDs in the original local
    extraction, while C/S contain the mapped cells ordered by those IDs.  The
    translation is only identifiable when every C/S column has one finite ID
    in the corresponding mapping column.  Refuse ambiguous files rather than
    silently selecting a neighbouring neuron.
    """
    mapping = np.asarray(mapping, dtype=float)
    rows = mapping_rows_for_columns(mapping, sub_index, n_columns)
    column_by_row = {int(row): column for column, row in enumerate(rows)}
    try:
        return np.asarray([column_by_row[int(row)] for row in global_ids], dtype=int)
    except KeyError as error:
        raise ValueError(
            f"la celula global {error.args[0]} no esta en la subsesion {sub_index + 1}"
        ) from error
