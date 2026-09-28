#!/usr/bin/env python3
"""Combine phasewise worker summaries without altering their original files."""

import argparse
import csv
import os
import time
from pathlib import Path


MAIN = Path('/mnt/NAS/Tomas/results/stoixeion/phasewise')
PARALLEL = Path('/mnt/NAS/Tomas/Miniscope-Data/results/stoixeion/active_parallel')
SOURCES = [
    MAIN / 'phase_summary.csv',
    PARALLEL / 'R004' / 'phase_summary.csv',
    PARALLEL / 'R005' / 'phase_summary.csv',
    PARALLEL / 'R006' / 'phase_summary.csv',
]
DESTINATION = MAIN / 'phase_summary_combined.csv'
KEY = ('animal', 'day_name', 'session', 'phase')


def combine():
    rows = {}
    fields = None
    for source in SOURCES:
        if not source.exists():
            continue
        with source.open(newline='') as handle:
            reader = csv.DictReader(handle)
            if fields is None:
                fields = reader.fieldnames
            for row in reader:
                rows.setdefault(tuple(row[name] for name in KEY), row)
    if fields is None:
        return 0
    temporary = DESTINATION.with_suffix('.csv.tmp')
    with temporary.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows[key] for key in sorted(rows))
    os.replace(temporary, DESTINATION)
    return len(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--watch-pid', type=int)
    args = parser.parse_args()
    while True:
        count = combine()
        print(f'{count} phase rows in {DESTINATION}', flush=True)
        if args.watch_pid is None or not Path(f'/proc/{args.watch_pid}').exists():
            break
        time.sleep(60)


if __name__ == '__main__':
    main()
