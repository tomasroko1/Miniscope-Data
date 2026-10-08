"""Object-zone enrichment with time-preserving circular-shift nulls."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from core.place_fields import gaussian_rate_maps_matlab, map_correlation
from .common import activity, bh_q, shift_choices, load_object_coords


def proximity_index(inside, outside):
    return np.divide(inside-outside, inside+outside, out=np.full_like(inside, np.nan),
                     where=(inside+outside)>0)


def zone_rates(phase, values, coords, radius, min_speed):
    xy = coords[['x_cm', 'y_cm']].to_numpy(float)
    zones = (phase['x'][:, None]-xy[None, :, 0])**2 + (phase['y'][:, None]-xy[None, :, 1])**2 <= radius**2
    valid = phase['valid'] & (phase['speed'] >= min_speed)
    inside = zones.any(axis=1) & valid
    outside = ~zones.any(axis=1) & valid
    occ_in = phase['dt'][inside].sum(); occ_out = phase['dt'][outside].sum()
    rates_in = values[inside].sum(axis=0)/occ_in if occ_in else np.full(values.shape[1], np.nan)
    rates_out = values[outside].sum(axis=0)/occ_out if occ_out else np.full(values.shape[1], np.nan)
    per = np.full((len(coords), values.shape[1]), np.nan)
    for i in range(len(coords)):
        mask = zones[:, i] & valid
        duration = phase['dt'][mask].sum()
        if duration:
            per[i] = values[mask].sum(axis=0)/duration
    return rates_in, rates_out, per, inside, outside, occ_in, occ_out


def analyse_phase(phase, coords, opts, rng):
    values = activity(phase['S'], opts.activity, opts.threshold_sd)
    values[~np.isfinite(values)] = 0
    rin, rout, per, inside, outside, oi, oo = zone_rates(phase, values, coords, opts.radius, opts.object_min_speed)
    observed = proximity_index(rin, rout)
    shifts = shift_choices(len(values), np.median(phase['dt']), rng, opts.object_shuffles)
    # FFT computes sums for all circular shifts at once for every cell.
    transformed = np.fft.rfft(values, axis=0)
    def shifted_sum(mask):
        return np.fft.irfft(np.conj(np.fft.rfft(mask.astype(float)))[:, None]*transformed,
                           n=len(values), axis=0)[shifts]
    null = proximity_index(shifted_sum(inside)/oi, shifted_sum(outside)/oo) if oi and oo else np.full((len(shifts), values.shape[1]), np.nan)
    nevents = values[phase['valid'] & (phase['speed'] >= opts.object_min_speed)].sum(axis=0)
    reasons = np.full(values.shape[1], '', dtype=object)
    reasons[nevents < opts.object_min_events] = 'actividad insuficiente'
    if oi < opts.min_object_occupancy or oo < opts.min_object_occupancy:
        reasons[:] = 'ocupación insuficiente'
    reasons[~np.isfinite(observed) & (reasons == '')] = 'índice indefinido'
    p = (1 + (null >= observed[None, :]).sum(axis=0))/(len(shifts)+1)
    p[(reasons != '') | (opts.object_shuffles == 0)] = np.nan
    q = bh_q(p)
    rows = []
    for i, gid in enumerate(phase['ids']):
        preferred = coords.object_id.iloc[int(np.nanargmax(per[:, i]))] if np.isfinite(per[:, i]).any() else ''
        row = dict(phase=phase['phase'], global_id=int(gid), id_scope='global' if phase['aligned'] else 'local',
                   local_column=int(phase['columns'][i]), n_events=float(nevents[i]), rate_in=rin[i], rate_out=rout[i],
                   object_proximity_index=observed[i], preferred_object=preferred, p_value=p[i], q_value=q[i],
                   is_object_cell=bool(p[i] < opts.alpha and observed[i] > 0),
                   is_object_cell_fdr=bool(q[i] < opts.alpha and observed[i] > 0),
                   sig_in_OF_control=False, coords_are_placeholder=bool(coords.is_placeholder.any()),
                   radius_cm=opts.radius, occupancy_in_s=oi, occupancy_out_s=oo, exclusion_reason=reasons[i],
                   activity_mode=opts.activity, min_speed=opts.object_min_speed)
        for j, oid in enumerate(coords.object_id):
            row['rate_object_'+str(oid)] = per[j, i]
        rows.append(row)
    return pd.DataFrame(rows), dict(values=values, null=null, coords=coords, phase=phase)


def cell_map(cache, index):
    stored = cache.setdefault('maps', {})
    if index in stored:
        return stored[index]
    phase = cache['phase']; good = phase['valid'] & cache['speed_mask']
    maps, edges, _ = gaussian_rate_maps_matlab(phase['x'][good], phase['y'][good], phase['t'][good],
                        cache['values'][good, index], bin_size=1, sigma=4)
    stored[index] = (maps[0], edges)
    return stored[index]


DISPLACEMENT_COLUMNS = [
    'global_id', 'preferred_sample_object', 'displacement', 'present_in_TEST',
    'old_x_cm', 'old_y_cm', 'new_x_cm', 'new_y_cm', 'coords_are_placeholder',
    'rate_TEST_old', 'rate_TEST_new', 'follow_index', 'sample_test_map_correlation',
    'peak_distance_old_cm', 'peak_distance_new_cm', 'significant_in_TEST',
]

def displacement(table, caches, opts):
    if 'SAMPLE' not in caches or 'TEST' not in caches:
        return pd.DataFrame(columns=DISPLACEMENT_COLUMNS)
    sc, tc = caches['SAMPLE'], caches['TEST']
    sp, tp = sc['phase'], tc['phase']
    if not sp['aligned'] or not tp['aligned']:
        print('[OBJETOS] Displacement omitido: mapping no identifica células entre fases.', flush=True)
        return pd.DataFrame(columns=DISPLACEMENT_COLUMNS)
    test_ids = {int(g): i for i, g in enumerate(tp['ids'])}
    rows = []
    chosen = table[(table.phase == 'SAMPLE') & table.is_object_cell]
    for _, cell in chosen.iterrows():
        gid = int(cell.global_id)
        old = sc['coords'][sc['coords'].object_id == cell.preferred_object].iloc[0]
        matches = tc['coords'][tc['coords'].match_id == cell.preferred_object]
        if len(matches) > 1:
            raise ValueError('match_id TEST duplicado: correspondencia no única')
        oldxy = np.array([old.x_cm, old.y_cm], dtype=float)
        newxy = matches[['x_cm', 'y_cm']].to_numpy(float)[0] if len(matches) else None
        status = 'removed' if newxy is None else ('stayed' if np.linalg.norm(newxy-oldxy) < 2 else 'moved')
        row = dict(global_id=gid, preferred_sample_object=cell.preferred_object,
                   displacement=status, present_in_TEST=gid in test_ids,
                   old_x_cm=oldxy[0], old_y_cm=oldxy[1],
                   new_x_cm=newxy[0] if newxy is not None else np.nan,
                   new_y_cm=newxy[1] if newxy is not None else np.nan,
                   coords_are_placeholder=bool(sc['coords'].is_placeholder.any() or tc['coords'].is_placeholder.any()),
                   rate_TEST_old=np.nan, rate_TEST_new=np.nan, follow_index=np.nan,
                   sample_test_map_correlation=np.nan, peak_distance_old_cm=np.nan,
                   peak_distance_new_cm=np.nan, significant_in_TEST=False)
        if gid in test_ids:
            ti = test_ids[gid]; si = int(np.flatnonzero(sp['ids'] == gid)[0])
            def near_rate(xy):
                valid = tp['valid'] & tc['speed_mask']
                mask = valid & ((tp['x']-xy[0])**2+(tp['y']-xy[1])**2 <= opts.radius**2)
                duration = tp['dt'][mask].sum()
                return tc['values'][mask, ti].sum()/duration if duration >= opts.min_object_occupancy else np.nan
            ro = near_rate(oldxy); rn = near_rate(newxy) if newxy is not None else np.nan
            sm, _ = cell_map(sc, si); tm, edges = cell_map(tc, ti)
            centres = (edges[:-1]+edges[1:])/2
            peak = np.array([centres[j] for j in np.unravel_index(np.nanargmax(tm), tm.shape)])
            row.update(rate_TEST_old=ro, rate_TEST_new=rn,
                       follow_index=(rn-ro)/(rn+ro) if status == 'moved' and rn+ro > 0 else np.nan,
                       sample_test_map_correlation=map_correlation(sm, tm),
                       peak_distance_old_cm=np.linalg.norm(peak-oldxy),
                       peak_distance_new_cm=np.linalg.norm(peak-newxy) if newxy is not None else np.nan,
                       significant_in_TEST=bool(table[(table.phase == 'TEST') & (table.global_id == gid)].is_object_cell.iloc[0]))
        rows.append(row)
    return pd.DataFrame(rows, columns=DISPLACEMENT_COLUMNS)


def example_figures(table, caches, folder, opts):
    ranked = table[table.phase == 'SAMPLE'].sort_values(['p_value', 'object_proximity_index'], ascending=[True, False])
    for _, row in ranked.head(opts.n_examples).iterrows():
        gid = int(row.global_id)
        names = [p for p in ['OF1', 'SAMPLE', 'TEST', 'OF2'] if p in caches]
        fig, axes = plt.subplots(3, len(names), figsize=(3.3*len(names), 9), squeeze=False)
        for j, name in enumerate(names):
            cache = caches[name]; phase = cache['phase']
            if not phase['aligned'] and name != 'SAMPLE':
                for ax in axes[:, j]: ax.set_visible(False)
                continue
            found = np.flatnonzero(phase['ids'] == gid)
            if not len(found):
                for ax in axes[:, j]: ax.set_visible(False)
                continue
            i = found[0]; valid = phase['valid'] & cache['speed_mask']
            x, y = phase['x'], phase['y']; active = valid & (cache['values'][:, i] > 0)
            axes[0, j].plot(x[valid], y[valid], color='.8', lw=.3)
            axes[0, j].scatter(x[active], y[active], s=4, c='red')
            rate, edges = cell_map(cache, i)
            im = axes[1, j].imshow(rate.T, origin='lower', extent=[edges[0], edges[-1]]*2)
            cb = fig.colorbar(im, ax=axes[1, j], fraction=.046, pad=.04)
            cb.set_label('events/s' if opts.activity == 'events' else 'S/s', fontsize=8)
            cb.ax.tick_params(labelsize=7)
            for ax in axes[:2, j]:
                for _, obj in cache['coords'].iterrows():
                    ax.add_patch(Circle((obj.x_cm, obj.y_cm), opts.radius, fill=False, color='orange'))
                ax.set(xlim=(-25,25), ylim=(-25,25), aspect='equal')
            axes[0, j].set_title(name+' | célula '+str(gid))
            null = cache['null'][:, i]; axes[2, j].hist(null[np.isfinite(null)], bins=25, color='.6')
            obs = table[(table.phase == name) & (table.global_id == gid)].object_proximity_index.iloc[0]
            axes[2, j].axvline(obs, color='red'); axes[2, j].set_xlabel('Índice proximidad (nulo)')
        fig.tight_layout(); fig.savefig(folder / ('example_%04d.png' % gid), dpi=110); plt.close(fig)


def analyse_session(phases, animal, day, opts, rng, folder):
    coords = {name: load_object_coords(opts.objects_csv, animal, day, name) for name in ['SAMPLE', 'TEST']}
    if not set(coords['TEST'].match_id) - {''} <= set(coords['SAMPLE'].object_id):
        raise ValueError('match_id TEST no existe en SAMPLE')
    nonempty = coords['TEST'].match_id[coords['TEST'].match_id != '']
    if nonempty.duplicated().any():
        raise ValueError('match_id TEST duplicado')
    tables, caches = [], {}
    for phase in phases:
        name = phase['phase']
        if name not in ['OF1', 'SAMPLE', 'TEST', 'OF2']:
            continue
        table, cache = analyse_phase(phase, coords['TEST' if name == 'TEST' else 'SAMPLE'], opts, rng)
        cache['speed_mask'] = phase['speed'] >= opts.object_min_speed
        tables.append(table); caches[name] = cache
    if not tables:
        raise ValueError('sesión sin fases de objetos reconocibles')
    table = pd.concat(tables, ignore_index=True)
    if all(p['aligned'] for p in phases):
        control_ids = set(table[table.phase.isin(['OF1', 'OF2']) & table.is_object_cell].global_id)
        table['sig_in_OF_control'] = table.global_id.isin(control_ids)
    else:
        table['sig_in_OF_control'] = pd.NA
    displaced = displacement(table, caches, opts)
    example_figures(table, caches, folder, opts)
    return table, displaced
