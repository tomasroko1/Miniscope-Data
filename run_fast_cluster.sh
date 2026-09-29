#!/usr/bin/env bash
# RUN_FAST_CLUSTER.SH
# Runs the fast Stoixeion analysis (99 shuffles, no figures, phasewise + global)
# for the 9 focal sessions: R004, R005, R006 in HabL, SD CNO, and SD VEH.

set -euo pipefail

REPO_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DATA_DIR=${1:-"/mnt/NAS/Miniscopes/Reg_CA1/DataBase"}
RESULTS_ROOT="/mnt/NAS/Tomas/results/stoixeion"
OUT_DIR=${2:-"${RESULTS_ROOT}/runs/$(date +%Y-%m-%d_%H%M%S)_fast_9_focal_sessions"}
DEFAULT_STOIXEION="${REPO_DIR}/external/Stoixeion"
if [[ -f "/mnt/NAS/Tomas/Stoixeion/Stoixeion.m" ]]; then
  DEFAULT_STOIXEION="/mnt/NAS/Tomas/Stoixeion"
fi
STOIXEION_DIR=${3:-"${DEFAULT_STOIXEION}"}
MATLAB_BIN=${MATLAB_BIN:-"/usr/local/MATLAB/R2017a/bin/matlab"}

if [[ -e "${OUT_DIR}" ]]; then
  echo "El directorio de salida ya existe: ${OUT_DIR}" >&2
  exit 2
fi
mkdir -p "${OUT_DIR}"
echo "=========================================================="
echo "Iniciando corrida rapida de Stoixeion (99 shuffles, sin PNGs)"
echo "Repo:       ${REPO_DIR}"
echo "Data:       ${DATA_DIR}"
echo "Stoixeion:  ${STOIXEION_DIR}"
echo "Destino:    ${OUT_DIR}"
echo "=========================================================="

# Ensure Stoixeion dependency is ready
if [[ ! -f "${STOIXEION_DIR}/Stoixeion.m" ]]; then
  echo "Preparando dependencia Stoixeion..."
  bash "${REPO_DIR}/setup_stoixeion.sh" "${STOIXEION_DIR}"
fi

ANIMALS=("R004" "R005" "R006")
PIDS=()

echo "Lanzando 3 workers en paralelo (uno por animal)..."
for ANIMAL in "${ANIMALS[@]}"; do
  WORKER_OUT="${OUT_DIR}/workers/${ANIMAL}"
  LOG_FILE="${OUT_DIR}/${ANIMAL}_run.log"
  mkdir -p "${WORKER_OUT}"
  
  echo "  -> Lanzando worker para ${ANIMAL}..."
  
  # Run worker in background (single-line -r prevents MATLAB syntax errors with newlines)
  (
    "${MATLAB_BIN}" -nodisplay -nodesktop -nosplash -r "try, cd('${REPO_DIR}'); addpath(pwd); run_stoixeion_fast('${ANIMAL}', '${DATA_DIR}', '${STOIXEION_DIR}', '${WORKER_OUT}', false); exit(0); catch ME, disp(getReport(ME, 'extended')); exit(1); end;" > "${LOG_FILE}" 2>&1
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

PYTHON_BIN=${PYTHON_BIN:-}
if [[ -z "${PYTHON_BIN}" ]]; then
  for CANDIDATE in /mnt/NAS/Tomas/results/minian_registration_env_20260929/bin/python python3 python; do
    if command -v "${CANDIDATE}" >/dev/null 2>&1 && \
       "${CANDIDATE}" -c 'import matplotlib, pandas' >/dev/null 2>&1; then
      PYTHON_BIN=${CANDIDATE}
      break
    fi
  done
fi

if [[ -n "${PYTHON_BIN}" ]]; then
  echo "Consolidando tablas y generando figuras de publicacion..."
  "${PYTHON_BIN}" "${REPO_DIR}/plot_publication_ensembles.py" \
      --results-dir "${OUT_DIR}" \
      --output-dir "${OUT_DIR}/figures" || echo "[AVISO] La generacion de graficos fallo. Las tablas estan intactas en ${OUT_DIR}/workers/"
else
  echo "[AVISO] No se encontro Python con pandas y matplotlib. Las tablas estan guardadas en ${OUT_DIR}/workers/"
fi

echo "=========================================================="
echo "Corrida completada exitosamente."
echo "Tablas consolidadas: ${OUT_DIR}/tables/"
echo "Figuras limpias:     ${OUT_DIR}/figures/"
echo "=========================================================="
if [[ "${OUT_DIR}" == "${RESULTS_ROOT}/runs/"* ]]; then
  ln -sfn "${OUT_DIR}" "${RESULTS_ROOT}/current"
fi
