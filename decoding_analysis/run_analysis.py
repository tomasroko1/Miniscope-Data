"""CLI: python -m decoding_analysis.run_analysis {decode,objects,all}."""
import argparse
import json
from pathlib import Path
import traceback
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
from core.data_loader import get_data_dir, list_sessions, load_session
from .common import session_phases, category
from .decoding import bin_phase, decode_phase
from .object_cells import analyse_session


def plot_nulls(nulls, rows, path):
    fig, axes = plt.subplots(1, len(rows), figsize=(3.4*len(rows), 3), squeeze=False)
    for ax, null, row in zip(axes[0], nulls, rows):
        ax.hist(null, bins=25, color='.6', label='shuffles')
        ax.axvline(row['median_error_cm'], color='red', label='observado')
        ax.axvline(row['chance_median_error_cm'], color='blue', ls=':', label='azar')
        ax.set(title=row['phase'], xlabel='Error mediano (cm)')
    axes[0,0].legend(fontsize=7); fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)


def decoding_reports(table, root):
    table.to_csv(root/'decoding_summary.csv', index=False)
    deltas = []
    for (animal, session), group in table.groupby(['animal', 'session_file']):
        errors = group.set_index('phase').median_error_cm
        z = group.set_index('phase').shuffle_z
        if all(p in errors.index for p in ['SAMPLE', 'TEST', 'OF1', 'OF2']):
            obj = errors[['SAMPLE','TEST']].mean(); of = errors[['OF1','OF2']].mean()
            deltas.append(dict(animal=animal, session_file=session, treatment=group.treatment.iloc[0],
                               task=group.task.iloc[0], objects_error_cm=obj, of_error_cm=of, delta_cm=obj-of,
                               objects_shuffle_z=z[['SAMPLE','TEST']].mean(),
                               of_shuffle_z=z[['OF1','OF2']].mean(),
                               delta_shuffle_z=z[['SAMPLE','TEST']].mean()-z[['OF1','OF2']].mean()))
    delta = pd.DataFrame(deltas, columns=['animal','session_file','treatment','task','objects_error_cm','of_error_cm','delta_cm','objects_shuffle_z','of_shuffle_z','delta_shuffle_z'])
    delta.to_csv(root/'decoding_objects_vs_of.csv', index=False)
    tests = []
    for name, values in [('overall', delta.delta_cm)] + [(t, g.delta_cm) for t,g in delta.groupby('treatment')]:
        p, statistic, note = np.nan, np.nan, 'n < 3'
        if len(values) >= 3:
            if np.all(values == 0): p, statistic = 1., 0.
            else:
                res = wilcoxon(values); p, statistic = res.pvalue, res.statistic
            note = 'Sesiones como pares; animales repetidos no son independientes'
        tests.append(dict(treatment=name, n_sessions=len(values), statistic=statistic, p_value=p, note=note))
    pd.DataFrame(tests).to_csv(root/'decoding_objects_vs_of_wilcoxon.csv', index=False)
    figures = root/'figures'; figures.mkdir(exist_ok=True)
    for label, tasks, phases in [('objects', ['SD','XsS'], ['OF1','SAMPLE','TEST','OF2']),
                                 ('HabL', ['HabL'], ['OF1','OF2','OF3','OF4']),
                                 ('other', ['other'], list(table[table.task == 'other'].phase.unique()))]:
        subset = table[table.task.isin(tasks)]
        if subset.empty: continue
        treatments = list(subset.treatment.unique())
        fig, axes = plt.subplots(1, len(treatments), figsize=(5*len(treatments),4), squeeze=False)
        for ax, treatment in zip(axes[0], treatments):
            for _, group in subset[subset.treatment == treatment].groupby(['animal','session_file']):
                group = group.set_index('phase').reindex(phases)
                ax.plot(range(len(phases)), group.median_error_cm, 'o-', alpha=.55, color='black')
                ax.plot(range(len(phases)), group.shuffle_median_of_medians_cm, '--', alpha=.4, color='orange')
                ax.plot(range(len(phases)), group.chance_median_error_cm, ':', alpha=.4, color='blue')
            ax.set(xticks=range(len(phases)), xticklabels=phases, ylabel='Error mediano (cm)', title=treatment)
        fig.suptitle('Negro: observado; naranja: shuffle; azul: azar')
        fig.tight_layout(); fig.savefig(figures/('phase_errors_'+label+'.png'), dpi=130); plt.close(fig)
    if not delta.empty:
        fig, ax = plt.subplots(figsize=(4,4))
        for i, (t, group) in enumerate(delta.groupby('treatment')):
            ax.scatter(np.full(len(group), i), group.delta_cm, label=t)
        names = sorted(delta.treatment.unique())
        ax.set(xticks=range(len(names)), xticklabels=names, ylabel='Objetos − OF: error (cm)')
        ax.axhline(0, color='.5'); fig.tight_layout(); fig.savefig(figures/'objects_vs_of.png', dpi=130); plt.close(fig)


def parser():
    p = argparse.ArgumentParser(description='Decodificación y células de objetos')
    p.add_argument('analysis', choices=['decode','objects','all'])
    p.add_argument('--data-dir', default=get_data_dir()); p.add_argument('--out-dir', default='results')
    p.add_argument('--animal'); p.add_argument('--session')
    p.add_argument('--objects-csv', default=str(Path(__file__).with_name('object_coords.csv')))
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--activity', choices=['events','S'], default='events')
    p.add_argument('--threshold-sd', type=float, default=3)
    p.add_argument('--time-bin', type=float, default=.5)
    p.add_argument('--decoder-bin', type=float, default=2)
    p.add_argument('--sigma', type=float, default=4)
    p.add_argument('--k-folds', type=int, default=5)
    p.add_argument('--fold-gap-s', type=float, default=2)
    p.add_argument('--prior', choices=['uniform','occupancy'], default='uniform')
    p.add_argument('--cells', choices=['common','all'], default='common')
    p.add_argument('--equalize-duration', action='store_true')
    p.add_argument('--save-predictions', action='store_true')
    p.add_argument('--min-speed', type=float, default=None, help='Default: decode=2, objetos=0 cm/s')
    p.add_argument('--min-events', type=float, default=None, help='Default: decode=1 en training, objetos=10')
    p.add_argument('--n-shuffles', type=int, default=None, help='Sobrescribe ambos defaults')
    p.add_argument('--decode-shuffles', type=int, default=100)
    p.add_argument('--object-shuffles', type=int, default=500)
    p.add_argument('--arena-tolerance-cm', type=float, default=2)
    p.add_argument('--radius', type=float, default=5)
    p.add_argument('--min-object-occupancy', type=float, default=2)
    p.add_argument('--alpha', type=float, default=.05)
    p.add_argument('--n-examples', type=int, default=6)
    p.add_argument('--verbose', action='store_true')
    return p


def main():
    p = parser(); opts = p.parse_args()
    opts.decode_min_speed = 2 if opts.min_speed is None else opts.min_speed
    opts.object_min_speed = 0 if opts.min_speed is None else opts.min_speed
    opts.decode_min_events = 1 if opts.min_events is None else opts.min_events
    opts.object_min_events = 10 if opts.min_events is None else opts.min_events
    if opts.n_shuffles is not None:
        opts.decode_shuffles = opts.object_shuffles = opts.n_shuffles
    if min(opts.time_bin, opts.decoder_bin, opts.sigma, opts.radius) <= 0 or opts.k_folds < 2:
        p.error('bins/sigma/radio > 0 y k-folds >= 2')
    if min(opts.decode_shuffles, opts.object_shuffles, opts.n_examples, opts.fold_gap_s,
           opts.decode_min_speed, opts.object_min_speed, opts.decode_min_events,
           opts.object_min_events, opts.min_object_occupancy, opts.arena_tolerance_cm) < 0 or not 0 < opts.alpha < 1:
        p.error('parámetros negativos o alpha inválida')
    if not np.isfinite(opts.arena_tolerance_cm):
        p.error('arena-tolerance-cm debe ser finita')
    root = Path(opts.out_dir)
    root.mkdir(parents=True, exist_ok=True)
    with (root/'run_config.json').open('w') as handle:
        json.dump(vars(opts), handle, indent=2, ensure_ascii=False)
        handle.write('\n')
    session_file = opts.session
    if session_file and not session_file.endswith('.mat'):
        session_file += '.mat'
    sessions = [s for s in list_sessions(opts.data_dir) if
                (not opts.animal or s['animal'] == opts.animal) and
                (not session_file or s['session_file'] == session_file)]
    if not sessions: raise SystemExit('No se encontraron sesiones')
    rng = np.random.default_rng(opts.seed)
    decoded, objects, displaced, failures = [], [], [], []
    for session in sessions:
        animal, filename = session['animal'], session['session_file']; stem = Path(filename).stem
        print('Procesando %s / %s' % (animal, filename), flush=True)
        try:
            data = load_session(animal, filename, opts.data_dir)
            day, phases = session_phases(data, arena_tolerance_cm=opts.arena_tolerance_cm)
            task, treatment = category(day)
            metadata = dict(animal=animal, session_file=filename, day_name=day, task=task, treatment=treatment)
        except Exception as err:
            print('[OMITIDA] %s' % err, flush=True); failures.append(dict(**session, analysis='load', reason=str(err))); continue
        for mode in ['decode', 'objects']:
            if opts.analysis not in [mode, 'all'] or (mode == 'objects' and task not in ['SD','XsS']): continue
            try:
                folder = root/('decoding' if mode == 'decode' else 'object_cells')/animal/stem
                folder.mkdir(parents=True, exist_ok=True)
                if mode == 'decode':
                    _, dp = session_phases(data, common=opts.cells == 'common', arena_tolerance_cm=opts.arena_tolerance_cm)
                    binned = [bin_phase(ph, opts.time_bin) for ph in dp]
                    eligible = [b['valid'] & (b['speed'] >= opts.decode_min_speed) for b in binned]
                    if opts.equalize_duration:
                        minimum = min(mask.sum() for mask in eligible)
                        for mask in eligible: mask[np.flatnonzero(mask)[minimum:]] = False
                    rows, nulls = [], []
                    for ph, b, mask in zip(dp, binned, eligible):
                        metrics, predictions, null = decode_phase(b, opts, rng, mask)
                        row = dict(metadata, frac_frames_projected=ph['frac_frames_projected'],
                                   frac_frames_outside_arena=ph['frac_frames_outside_arena'], phase=ph['phase'], has_objects=ph['phase'] in ['SAMPLE','TEST'], **metrics)
                        rows.append(row); nulls.append(null)
                        if opts.save_predictions:
                            predictions.to_csv(folder/('predictions_'+ph['phase'].replace('/','_')+'.csv'), index=False)
                        pd.DataFrame({'median_error_cm':null}).to_csv(folder/('shuffle_'+ph['phase'].replace('/','_')+'.csv'), index=False)
                        print('  %s: %.2f cm; nulo p5 %.2f; azar %.2f' %
                              (ph['phase'], metrics['median_error_cm'], metrics['shuffle_p5_cm'], metrics['chance_median_error_cm']), flush=True)
                    pd.DataFrame(rows).to_csv(folder/('decoding_'+stem+'.csv'), index=False)
                    figures = root/'decoding'/'figures'; figures.mkdir(exist_ok=True)
                    plot_nulls(nulls, rows, figures/(animal+'_'+stem+'_null.png'))
                    decoded.extend(rows)
                else:
                    table, disp = analyse_session(phases, animal, day, opts, rng, folder)
                    for k, v in metadata.items(): table[k] = v; disp[k] = v
                    table.to_csv(folder/('object_cells_'+stem+'.csv'), index=False)
                    disp.to_csv(folder/('object_displacement_'+stem+'.csv'), index=False)
                    objects.append(table); displaced.append(disp)
                    print('  SAMPLE: %d células significativas' % table[table.phase == 'SAMPLE'].is_object_cell.sum(), flush=True)
            except Exception as err:
                print('[OMITIDA %s] %s/%s: %s' % (mode, animal, filename, err), flush=True)
                failures.append(dict(**session, analysis=mode, reason=str(err)))
                if opts.verbose: traceback.print_exc()
    root.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(failures, columns=['animal','session_file','analysis','reason']).to_csv(root/'analysis_failures.csv', index=False)
    if decoded: decoding_reports(pd.DataFrame(decoded), root/'decoding')
    if objects:
        table = pd.concat(objects, ignore_index=True)
        table.to_csv(root/'object_cells'/'object_cells_summary.csv', index=False)
        disp = pd.concat(displaced, ignore_index=True)
        disp.to_csv(root/'object_cells'/'object_displacement_summary.csv', index=False)
        counts = table.groupby(['animal','session_file','task','treatment','phase']).agg(
            n_cells=('global_id','size'), n_eligible=('p_value','count'),
            n_object_cells=('is_object_cell','sum'), n_object_cells_strict=('is_object_cell_strict','sum'), n_object_cells_fdr=('is_object_cell_fdr','sum')).reset_index()
        counts['fraction_object_cells_strict'] = counts.n_object_cells_strict/counts.n_cells
        counts['fraction_object_cells'] = counts.n_object_cells/counts.n_cells
        counts.to_csv(root/'object_cells'/'object_cell_counts.csv', index=False)
        for (a, f), group in counts.groupby(['animal','session_file']):
            group.to_csv(root/'object_cells'/a/Path(f).stem/'object_cell_counts.csv', index=False)
    print('Listo. Fallos registrados: %d' % len(failures), flush=True)
    if failures: raise SystemExit(1)


if __name__ == '__main__': main()
