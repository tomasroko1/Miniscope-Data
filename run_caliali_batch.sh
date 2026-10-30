#!/usr/bin/env bash
set -e

# Ir a la raíz del repositorio
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_DIR}"

# Seleccionar el Python del entorno virtual (.venv) que tiene matplotlib, scipy, etc.
if [ -f "${REPO_DIR}/.venv/bin/python" ]; then
    PYTHON="${REPO_DIR}/.venv/bin/python"
elif [ -f "${REPO_DIR}/.venv/bin/python3" ]; then
    PYTHON="${REPO_DIR}/.venv/bin/python3"
elif command -v python3 >/dev/null 2>&1; then
    PYTHON="python3"
else
    PYTHON="python"
fi

echo "=========================================================="
echo "Usando intérprete: ${PYTHON}"
${PYTHON} -c "import sys; print('Python version: ' + sys.version.split()[0])"
${PYTHON} -c "import matplotlib; print('Matplotlib version: ' + matplotlib.__version__)" 2>/dev/null || echo "[AVISO] matplotlib no encontrado en este python"
echo "=========================================================="

DATA_DIR="/mnt/NAS/Miniscopes/Reg_CA1/DataBase"
CALIALI_DIR="/mnt/NAS/Miniscopes/Reg_CA1/GDi"
OUT_DIR="${REPO_DIR}/results/Mapas_por_tipo_caliali"

echo "Data (tracking):   ${DATA_DIR}"
echo "CaliAli (senales): ${CALIALI_DIR}"
echo "Destino salida:    ${OUT_DIR}"
echo "=========================================================="

mkdir -p "${OUT_DIR}"

# 1. Batch principal de mapas espaciales
echo ">> Ejecutando run_batch.py..."
${PYTHON} run_batch.py \
    --data-dir "${DATA_DIR}" \
    --caliali-dir "${CALIALI_DIR}" \
    --out-dir "${OUT_DIR}" "$@"

# 2. Diagnóstico visual para sesiones excluidas (R004-XsS-CNO y R005-XsS-VEH)
echo ">> Ejecutando diagnostico de sesiones excluidas..."
${PYTHON} generar_diagnostico_dlc_excluidos.py \
    --data-dir "${DATA_DIR}" \
    --out-dir "${OUT_DIR}"

echo "=========================================================="
echo "Listo! Proceso finalizado correctamente."
echo "Resultados en: ${OUT_DIR}"
echo "=========================================================="
