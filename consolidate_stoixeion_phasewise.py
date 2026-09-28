#!/usr/bin/env python3
"""Keep one phasewise result tree while the R004 recovery worker finishes."""

import argparse
import csv
import os
import shutil
import time
from pathlib import Path


ROOT = Path('/mnt/NAS/Tomas/results/stoixeion/phasewise')
SOURCES = ROOT / 'sources'
R004_WORKER = Path('/mnt/NAS/Tomas/Miniscope-Data/results/stoixeion/active_parallel/R004')
NAMES = ('phase_summary', 'core_members', 'core_overlap', 'core_shuffle_summary',
         'ensemble_activity', 'ensemble_members', 'event_thresholds', 'singular_values')
PHASE_KEY = ('animal', 'day_name', 'session', 'phase')


def source_dirs():
    return [SOURCES / 'R004_original',
            SOURCES / 'R004_resume' if (SOURCES / 'R004_resume').exists() else R004_WORKER,
            SOURCES / 'R005', SOURCES / 'R006']


def read_rows(path):
    if not path.exists():
        return None, []
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


def consolidate():
    sources = source_dirs()
    preferred = {}
    for source in sources:
        _, rows = read_rows(source / 'phase_summary.csv')
        for row in rows:
            preferred.setdefault(tuple(row[key] for key in PHASE_KEY), source)

    for name in NAMES:
        fields = None
        merged = []
        seen = set()
        for source in sources:
            current_fields, rows = read_rows(source / f'{name}.csv')
            if fields is None and current_fields:
                fields = current_fields
            for row in rows:
                if 'phase' in row and preferred.get(tuple(row[key] for key in PHASE_KEY)) != source:
                    continue
                identity = tuple(row.get(key, '') for key in current_fields)
                if identity in seen:
                    continue
                seen.add(identity)
                merged.append(row)
        if fields is None:
            continue
        output = ROOT / f'{name}.csv'
        temporary = output.with_suffix('.csv.tmp')
        with temporary.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(merged)
        os.replace(temporary, output)
    return len(preferred)


def finalize_r004():
    destination = SOURCES / 'R004_resume'
    destination.mkdir(parents=True, exist_ok=True)
    for item in R004_WORKER.glob('*.csv'):
        shutil.move(str(item), str(destination / item.name))
    protocol = R004_WORKER / 'protocol.txt'
    if protocol.exists():
        shutil.move(str(protocol), str(destination / protocol.name))

    worker_figures = R004_WORKER / 'R004'
    if worker_figures.exists():
        for source in worker_figures.rglob('*'):
            if not source.is_file():
                continue
            relative = source.relative_to(worker_figures)
            target = ROOT / 'R004' / relative
            if target.exists():
                target = destination / 'duplicate_figures' / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(source), str(target))
    consolidate()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--watch-pid', type=int)
    parser.add_argument('--finalize', action='store_true',
                        help='Move completed R004 worker files into the canonical folder now')
    args = parser.parse_args()
    if args.finalize:
        finalize_r004()
        print('R004 recovery files moved into phasewise', flush=True)
        return
    while True:
        count = consolidate()
        print(f'{count} phases consolidated in {ROOT}', flush=True)
        if args.watch_pid is None:
            break
        if not Path(f'/proc/{args.watch_pid}').exists():
            finalize_r004()
            print('R004 recovery files moved into phasewise', flush=True)
            break
        time.sleep(60)


if __name__ == '__main__':
    main()
