#!/usr/bin/env bash
# RUN_FAST_CLUSTER.SH
# Runs the fast Stoixeion analysis (99 shuffles, no figures, phasewise + global)
# for the 9 focal sessions: R004, R005, R006 in HabL, SD CNO, and SD VEH.

set -euo pipefail

REPO_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DATA_DIR=${1:-"/mnt/NAS/Miniscopes/Reg_CA1/DataBase"}
STOIXEION_DIR=${2:-"${REPO_DIR}/external/Stoixeion"}
OUT_DIR=${3:-"/mnt/NAS/Tomas/results/stoixeion/fast_target9"}
MATLAB_BIN=${MATLAB_BIN:-"/usr/local/MATLAB/R2017a/bin/matlab"}

mkdir -p "${OUT_DIR}"
echo "=========================================================="
echo "Iniciando corrida rapida de Stoixeion (99 shuffles, sin PNGs)"
echo "Repo:       ${REPO_DIR}"
echo "Data:       ${DATA_DIR}"
echo "Stoixeion:  ${STOIXEION_DIR}"
echo "Salida:     ${OUT_DIR}"
echo "=========================================================="

# Ensure Stoixeion dependency is ready
if [[ ! -f "${STOIXEION_DIR}/Stoixeion.m" ]]; then
  echo "Preparando dependencia Stoixeion..."
  bash "${REPO_DIR}/setup_stoixeion.sh" "${STOIXEION_DIR}"
fi

# The 9 target sessions:
# R004: 2026_06_18_merged (HabL), 2026_06_19_merged (T1_SD_CNO), 2026_06_21_merged (T2_SD_VEH_)
# R005: 2026_06_18_merged (HabL), 2026_06_19_merged (T1_SD_VEH), 2026_06_21_merged (T2_SD_CNO)
# R006: 2026_06_18_merged (HabL), 2026_06_20_merged (T1_SD_VEH), 2026_06_22_merged (T2_SD_CNO)

ANIMALS=("R004" "R005" "R006")
PIDS=()

echo "Lanzando 3 workers en paralelo (uno por animal)..."
for ANIMAL in "${ANIMALS[@]}"; do
  WORKER_OUT="${OUT_DIR}/${ANIMAL}_worker"
  LOG_FILE="${OUT_DIR}/${ANIMAL}_run.log"
  mkdir -p "${WORKER_OUT}"
  
  echo "  -> Lanzando worker para ${ANIMAL} (Log: ${LOG_FILE})..."
  
  # Run worker in background
  (
    "${MATLAB_BIN}" -nodisplay -nodesktop -nosplash -r "
      try
        addpath('${REPO_DIR}');
        run_stoixeion_fast('target9', '${DATA_DIR}', '${STOIXEION_DIR}', '${WORKER_OUT}', false);
        exit(0);
      catch ME
        disp(getReport(ME, 'extended'));
        exit(1);
      end;
    " > "${LOG_FILE}" 2>&1
  ) &
  PIDS+=($!)
done

echo "Esperando que completen los workers..."
FAILED=0
for i in "${!PIDS[@]}"; do
  PID="${PIDS[$i]}"
  ANIMAL="${ANIMALS[$i]}"
  if wait "${PID}"; then
    echo "  [OK] Worker ${ANIMAL} finalizo correctamente."
  else
    echo "  [ERROR] Worker ${ANIMAL} fallo. Revisar log: ${OUT_DIR}/${ANIMAL}_run.log"
    FAILED=$((FAILED + 1))
  fi
done

if [[ ${FAILED} -gt 0 ]]; then
  echo "Error: ${FAILED} workers fallaron." >&2
  exit 1
fi

echo "Consolidando tablas de salida..."
python3 "${REPO_DIR}/consolidate_stoixeion_phasewise.py" "${OUT_DIR}" || true

echo "=========================================================="
echo "Corrida completada exitosamente."
echo "Resultados consolidados en: ${OUT_DIR}"
echo "=========================================================="
