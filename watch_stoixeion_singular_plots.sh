#!/usr/bin/env bash
set -u

repo=/mnt/NAS/Tomas/Miniscope-Data
r004=/mnt/NAS/Tomas/results/stoixeion/phasewise
r005=${repo}/results/stoixeion/active_parallel/R005
r006=${repo}/results/stoixeion/active_parallel/R006
last_plot=0

while true; do
  latest=0
  for root in "$r004" "$r005" "$r006"; do
    if [[ -f "$root/singular_values.csv" ]]; then
      modified=$(stat -c %Y "$root/singular_values.csv")
      if (( modified > latest )); then latest=$modified; fi
    fi
  done

  if (( latest > last_plot )); then
    sleep 5
    /usr/local/MATLAB/R2017a/bin/matlab -nodisplay -nodesktop -nosplash -r \
      "try, cd('${repo}'); addpath(pwd); plot_stoixeion_singular_values('${r004}','${r005}','${r006}'); catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"
    if [[ $? -eq 0 ]]; then last_plot=$latest; fi
  fi

  running=false
  for pid in 57353 51929 52019; do
    if ps -p "$pid" -o args= 2>/dev/null | grep -q 'run_stoixeion_miniscope'; then
      running=true
      break
    fi
  done
  if [[ "$running" == false ]]; then break; fi
  sleep 60
done
