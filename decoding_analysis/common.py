"""Adaptación de sesiones, coordenadas y utilidades compartidas."""
import re
import warnings
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from core.data_loader import get_common_neurons, mapping_columns
from core.place_fields import threshold_event_frames


def category(day):
    # Mirrors run_batch._session_category without importing its plotting runner.
    m = re.fullmatch(r'T\d+_(SD|XsS)_(CNO|VEH)_?', day)
    return m.groups() if m else ('HabL' if day == 'HabL' else 'other', 'NA')


def phase_names(day, n):
    if n == 1:
        return [day or 'Registro']
    if day == 'HabL' and n == 4:
        return ['OF1', 'OF2', 'OF3', 'OF4']
    if ('_SD_' in day or '_XsS_' in day) and n == 4:
        return ['OF1', 'SAMPLE', 'TEST', 'OF2']
    return ['Sub %d' % (i + 1) for i in range(n)]


def session_phases(data, common=False):
    subs = [data] if data['type'] == 'simple' else data['subsessions']
    day = str(getattr(data.get('sess'), 'day_name', ''))
    mapping = data.get('mapping')
    aligned = False
    ids, cols = [], []
    if mapping is not None:
        try:
            mapping = np.asarray(mapping, dtype=float)
            if mapping.ndim != 2 or mapping.shape[1] != len(subs):
                raise ValueError('forma de mapping incompatible')
            shared = get_common_neurons(mapping)
            for i, sub in enumerate(subs):
                present = np.flatnonzero(np.isfinite(mapping[:, i]))
                # Validate every local column even if only common cells requested.
                mapping_columns(mapping, present, i, sub['S'].shape[1])
                chosen = shared if common else present
                ids.append(chosen)
                cols.append(mapping_columns(mapping, chosen, i, sub['S'].shape[1]))
            if common and len(shared) == 0:
                raise ValueError('no hay células comunes')
            aligned = True
        except ValueError as err:
            print('[MAPPING] %s; células locales por fase, sin alineación.' % err, flush=True)
            ids, cols = [], []
    if not ids:
        if data['type'] == 'merged' and mapping is None:
            print('[MAPPING] mapping ausente; células locales por fase.', flush=True)
        for sub in subs:
            cols.append(np.arange(sub['S'].shape[1]))
            ids.append(cols[-1].copy())
    result = []
    for name, sub, gids, columns in zip(phase_names(day, len(subs)), subs, ids, cols):
        x, y, t = [np.asarray(sub[k], dtype=float).ravel() for k in ('x', 'y', 't')]
        S = np.asarray(sub['S'], dtype=float)[:, columns]
        if len(t) < 2 or np.any(~np.isfinite(t)) or np.any(np.diff(t) <= 0):
            raise ValueError('timestamps inválidos/no crecientes')
        dt = np.r_[np.diff(t), np.median(np.diff(t))]
        speed = sub.get('v')
        if speed is None:
            speed = np.hypot(np.gradient(x, t), np.gradient(y, t))
            speed = gaussian_filter1d(speed, max(1, .25 / np.median(dt)))
        else:
            speed = np.asarray(speed, dtype=float).ravel()
        valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(speed) & np.isfinite(S).all(axis=1)
        valid &= x*x + y*y <= 25**2 + 1e-6
        result.append(dict(phase=name, x=x, y=y, t=t, dt=dt, speed=speed,
                           S=S, valid=valid, ids=gids, columns=columns, aligned=aligned))
    return day, result


def activity(S, mode='events', threshold_sd=3, source=None):
    if mode == 'S':
        return np.maximum(S, 0)
    src = S if source is None else source
    return np.column_stack([threshold_event_frames(S[:, i], threshold_source=src[:, i],
                             threshold_sd=threshold_sd)[0] for i in range(S.shape[1])])


def shift_choices(n, dt, rng, count):
    if count == 0:
        return np.array([], dtype=int)
    minimum = int(np.ceil(30 / dt))
    if n < 2 * minimum + 1:
        raise ValueError('duración insuficiente para shifts de al menos 30 s en ambos sentidos')
    return rng.integers(minimum, n - minimum + 1, size=count)


def load_object_coords(csv_path, animal, day_name, phase):
    """Select one complete row set; never mix specificity levels within a phase."""
    table = pd.read_csv(csv_path, comment='#', keep_default_na=False)
    required = ['session_type', 'day_name', 'animal', 'phase', 'object_id', 'x_cm',
                'y_cm', 'match_id', 'is_placeholder', 'note']
    if any(k not in table for k in required):
        raise ValueError('CSV de coordenadas: faltan columnas requeridas')
    task = category(day_name)[0]
    table = table[(table.session_type == task) & (table.phase == phase) &
                  table.animal.isin(['', animal]) & table.day_name.isin(['', day_name])].copy()
    if table.empty:
        raise ValueError('sin coordenadas para %s/%s/%s' % (animal, day_name, phase))
    table['specificity'] = (table.day_name != '').astype(int) * 2 + (table.animal != '').astype(int)
    table = table[table.specificity == table.specificity.max()].copy()
    if table.object_id.duplicated().any():
        raise ValueError('object_id duplicado en el conjunto de coordenadas')
    for key in ['x_cm', 'y_cm', 'is_placeholder']:
        table[key] = pd.to_numeric(table[key], errors='raise')
    if not np.isfinite(table[['x_cm', 'y_cm']].to_numpy()).all():
        raise ValueError('coordenadas no finitas')
    if (table.x_cm**2 + table.y_cm**2 > 625).any():
        raise ValueError('objeto fuera de la arena')
    if not table.is_placeholder.isin([0, 1]).all():
        raise ValueError('is_placeholder debe ser 0/1')
    if table.is_placeholder.any():
        warnings.warn('¡¡¡ COORDENADAS PLACEHOLDER !!! %s / %s / %s: reemplazar por mediciones reales.' %
                      (animal, day_name, phase), RuntimeWarning)
    return table.drop(columns='specificity').reset_index(drop=True)


def bh_q(p):
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    good = np.flatnonzero(np.isfinite(p))
    order = good[np.argsort(p[good])]
    if len(order):
        q[order] = np.minimum(1, np.minimum.accumulate((p[order] * len(order) /
                              np.arange(1, len(order)+1))[::-1])[::-1])
    return q
