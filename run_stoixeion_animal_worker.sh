#!/usr/bin/env bash
set -euo pipefail

animal=${1:?Usage: run_stoixeion_animal_worker.sh R004|R005|R006}
case "$animal" in
  R004|R005|R006) ;;
  *) echo "Supported animals: R004, R005, R006" >&2; exit 2 ;;
esac

repo_dir=/mnt/NAS/Tomas/Miniscope-Data
source_dir=/mnt/NAS/Miniscopes/Reg_CA1/DataBase
stoixeion_dir=/mnt/NAS/Tomas/Stoixeion
input_dir="/tmp/stoixeion-phasewise-${animal}"
output_dir="${repo_dir}/results/stoixeion/active_parallel/${animal}"
mkdir -p "${input_dir}/${animal}" "${output_dir}"
if [[ "$animal" == R004 ]]; then
  # The original worker completed HabL and T1, then stopped during T2.
  ln -sfn "${source_dir}/${animal}/2026_06_21_merged.mat" \
    "${source_dir}/${animal}/2026_06_24_merged.mat" \
    "${source_dir}/${animal}/2026_06_30_merged.mat" \
    "${input_dir}/${animal}/"
else
  ln -sfn "${source_dir}/${animal}/"*.mat "${input_dir}/${animal}/"
fi

exec /usr/local/MATLAB/R2017a/bin/matlab -nodisplay -nodesktop -nosplash -r \
  "try, cd('${repo_dir}'); addpath(pwd); run_stoixeion_miniscope('phasewise','${input_dir}','${stoixeion_dir}','${output_dir}'); catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"
