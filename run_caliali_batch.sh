#!/usr/bin/env bash
set -e

# Ir a la raíz del repositorio
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_DIR}"

# Buscar un intérprete de Python que tenga matplotlib instalado
PYTHON=""
CANDIDATES=(
    "${REPO_DIR}/.venv/bin/python"
    "${REPO_DIR}/.venv/bin/python3"
    "/mnt/NAS/Tomas/results/minian_registration_env_20260929/bin/python"
    "/mnt/NAS/Tomas/.venv/bin/python"
    "${HOME}/.venv/bin/python"
    "${HOME}/miniconda3/bin/python"
    "${HOME}/anaconda3/bin/python"
    "python3"
    "python"
)

for CAND in "${CANDIDATES[@]}"; do
    if [ -x "${CAND}" ] || command -v "${CAND}" >/dev/null 2>&1; then
        if "${CAND}" -c "import matplotlib, scipy, h5py" >/dev/null 2>&1; then
            PYTHON="${CAND}"
            break
        elif "${CAND}" -c "import matplotlib" >/dev/null 2>&1; then
            PYTHON="${CAND}"
            break
        fi
    fi
done

# Si ningún Python tiene matplotlib, instalarlo en el espacio de usuario
if [ -z "${PYTHON}" ]; then
    echo "Ningún Python del sistema tiene matplotlib instalado."
    echo "Instalando dependencias mínimas en el espacio de usuario (~/.local)..."
    python3 -m pip install --user matplotlib h5py scipy numpy || pip3 install --user matplotlib h5py scipy numpy || true
    PYTHON="python3"
fi

echo "=========================================================="
echo "Usando intérprete: ${PYTHON}"
${PYTHON} -c "import sys; print('Python: ' + sys.version.split()[0] + ' (' + sys.executable + ')')"
${PYTHON} -c "import matplotlib; print('Matplotlib: ' + matplotlib.__version__)"
${PYTHON} -c "import h5py; print('h5py: ' + h5py.__version__)"
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
