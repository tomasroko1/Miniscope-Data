#!/usr/bin/env python3
"""Replace Stoixeion's crowded activity plot from saved MAT and CSV data."""

import argparse
import csv
import os
import time
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


DATA = Path('/mnt/NAS/Miniscopes/Reg_CA1/DataBase')
CANONICAL = Path('/mnt/NAS/Tomas/results/stoixeion/phasewise')
WORKER = Path(__file__).resolve().parent / 'results/stoixeion/active_parallel/R004'


def rows(path):
    if not path.exists():
        return []
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


def phase_name(day, index):
    return ('OF1', 'OF2', 'OF3', 'OF4')[index] if day == 'HabL' else \
        ('OF1', 'SAMPLE', 'TEST', 'OF2')[index]


def plot_session(source, figure_root, animal, day, session, phase_rows, thresholds, events):
    mat_path = DATA / animal / session
    act = loadmat(str(mat_path), squeeze_me=True, struct_as_record=False,
                  variable_names=['act'])['act']
    phases = np.ravel(act.S)
    mapping = np.asarray(act.mapping, dtype=float) if hasattr(act, 'mapping') else None
    for phase_row in phase_rows:
        phase = phase_row['phase']
        output_dir = figure_root / animal / day / Path(session).stem / phase
        old_plot = output_dir / 'stoixeion_04.png'
        if not old_plot.exists():
            continue
        phase_index = next((i for i in range(4) if phase_name(day, i) == phase), None)
        if phase_index is None:
            continue
        key = (animal, day, session, phase)
        threshold_rows = sorted(thresholds.get(key, []),
                                key=lambda row: int(row['local_cell_index']))
        if not threshold_rows:
            print(f'Sin umbrales: {output_dir}; se retira la figura ilegible', flush=True)
            old_plot.unlink()
            continue
        signal = np.asarray(phases[phase_index])
        local_column = None
        if mapping is not None:
            valid_rows = np.flatnonzero(np.isfinite(mapping[:, phase_index]))
            sorted_rows = valid_rows[np.argsort(mapping[valid_rows, phase_index])]
            local_column = {int(global_row): index for index, global_row in enumerate(sorted_rows)}
        columns = []
        limits = []
        for row in threshold_rows:
            global_id = row['mapped_cell_id']
            if local_column is not None and global_id and global_id.lower() != 'nan':
                column = local_column[int(float(global_id))]
            else:
                column = int(row['local_cell_index']) - 1
            columns.append(column)
            limits.append(float(row['S_threshold']))
        if any(column < 0 or column >= signal.shape[1] for column in columns):
            raise ValueError(f'Columnas fuera de rango: {output_dir}')
        selected = signal[:, columns]
        population = np.count_nonzero(np.isfinite(selected) &
                                      (selected > np.asarray(limits)), axis=1)
        n_frames = len(population)
        n_factors = int(float(phase_row['n_ensembles']))
        fig, (ax_count, ax_raster) = plt.subplots(2, 1, sharex=True,
                                                   figsize=(13, 6.5),
                                                   gridspec_kw={'height_ratios': [2, 1]},
                                                   layout='constrained')
        fig.patch.set_facecolor('white')
        ax_count.plot(np.arange(1, n_frames + 1), population, color='black',
                      linewidth=0.45)
        ax_count.set_title('Actividad poblacional por frame')
        ax_count.set_ylabel('Neuronas activas')
        ax_count.set_ylim(bottom=0)
        for factor in range(1, n_factors + 1):
            frames = events.get(key, {}).get(factor, [])
            if frames:
                ax_raster.plot(frames, np.full(len(frames), factor), '|',
                               color='black', markersize=3.5, markeredgewidth=0.6,
                               linestyle='None')
        ax_raster.set_title('Apariciones de ensambles')
        ax_raster.set_xlabel('Frame')
        ax_raster.set_ylabel('Ensamble')
        ax_raster.set_xlim(1, max(n_frames, 2))
        ax_raster.set_ylim(0.5, max(n_factors + 0.5, 1.5))
        if n_factors:
            ax_raster.set_yticks(np.arange(1, n_factors + 1))
        else:
            ax_raster.set_yticks([])
            ax_raster.text(0.5, 0.5, 'Sin ensambles detectados',
                           transform=ax_raster.transAxes, ha='center', va='center')
        for axis in (ax_count, ax_raster):
            axis.spines[['top', 'right']].set_visible(False)
        temporary = output_dir / 'stoixeion_04.tmp.png'
        fig.savefig(temporary, dpi=160)
        plt.close(fig)
        os.replace(temporary, old_plot)
        (output_dir / 'stoixeion_05.png').unlink(missing_ok=True)
        print(f'Actualizada: {old_plot}', flush=True)


def process_source(source, figure_root, done):
    summaries = rows(source / 'phase_summary.csv')
    thresholds = defaultdict(list)
    for row in rows(source / 'event_thresholds.csv'):
        thresholds[tuple(row[k] for k in ('animal', 'day_name', 'session', 'phase'))].append(row)
    events = defaultdict(lambda: defaultdict(list))
    for row in rows(source / 'ensemble_activity.csv'):
        key = tuple(row[k] for k in ('animal', 'day_name', 'session', 'phase'))
        events[key][int(row['ensemble'])].append(int(row['frame']))
    sessions = defaultdict(list)
    for row in summaries:
        if row['phase'] == 'GLOBAL_CONCAT' or row['status'] != 'ok':
            continue
        key = tuple(row[k] for k in ('animal', 'day_name', 'session', 'phase'))
        if key in done:
            continue
        output = figure_root / key[0] / key[1] / Path(key[2]).stem / key[3] / 'stoixeion_04.png'
        if output.exists():
            sessions[key[:3]].append(row)
    for (animal, day, session), phase_rows in sessions.items():
        plot_session(source, figure_root, animal, day, session, phase_rows, thresholds, events)
        done.update((animal, day, session, row['phase']) for row in phase_rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--watch-pid', type=int)
    args = parser.parse_args()
    done = set()
    while True:
        process_source(CANONICAL, CANONICAL, done)
        process_source(WORKER, WORKER, done)
        if args.watch_pid is None or not Path(f'/proc/{args.watch_pid}').exists():
            break
        time.sleep(60)


if __name__ == '__main__':
    main()
