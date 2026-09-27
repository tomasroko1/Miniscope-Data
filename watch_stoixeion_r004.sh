#!/usr/bin/env bash
set -euo pipefail

original_pid=${1:?Usage: watch_stoixeion_r004.sh MATLAB_PID}
summary=/mnt/NAS/Tomas/results/stoixeion/phasewise/phase_summary.csv

while kill -0 "$original_pid" 2>/dev/null; do
  if [[ -f "$summary" ]] && awk -F, 'NR > 1 && $1 == "R005" {found=1; exit} END {exit !found}' "$summary"; then
    printf 'R004 finished; stopping original MATLAB process %s after it saved the first R005 phase.\n' "$original_pid"
    kill -TERM "$original_pid"
    exit 0
  fi
  sleep 15
done

printf 'Original MATLAB process %s exited before R005 appeared.\n' "$original_pid"
