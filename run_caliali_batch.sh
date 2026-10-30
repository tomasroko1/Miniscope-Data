#!/usr/bin/env bash
set -e

# Ir a la raíz del repositorio
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_DIR}"

# 1. Cargar hooks de conda si existen en el sistema
for CONDA_SH in \
    "${HOME}/miniconda3/etc/profile.d/conda.sh" \
    "${HOME}/anaconda3/etc/profile.d/conda.sh" \
    "/opt/conda/etc/profile.d/conda.sh" \
    "/usr/local/miniconda3/etc/profile.d/conda.sh"; do
    if [ -f "${CONDA_SH}" ]; then
        # shellcheck disable=SC1090
        source "${CONDA_SH}" 2>/dev/null || true
        break
    fi
done

PYTHON=""

# Si el usuario pasó un argumento que es un ejecutable de python o nombre de env
if [ -n "$1" ]; then
    if [ -x "$1" ]; then
        PYTHON="$1"
        shift
    elif command -v conda >/dev/null 2>&1 && conda env list | grep -q "^$1 "; then
        ENV_PATH=$(conda env list | grep "^$1 " | awk '{print $NF}')
        if [ -x "${ENV_PATH}/bin/python" ]; then
            PYTHON="${ENV_PATH}/bin/python"
            shift
        fi
    fi
fi

# 2. Si no se especificó directamente, recolectar candidatos
if [ -z "${PYTHON}" ]; then
    CANDIDATES=()

    # Si hay variable de entorno PYTHON_BIN
    if [ -n "${PYTHON_BIN}" ] && [ -x "${PYTHON_BIN}" ]; then
        CANDIDATES+=("${PYTHON_BIN}")
    fi

    # Si hay un entorno de Conda activado en este momento
    if [ -n "${CONDA_PREFIX}" ] && [ -x "${CONDA_PREFIX}/bin/python" ]; then
        CANDIDATES+=("${CONDA_PREFIX}/bin/python")
    fi

    # Registrar entornos de Conda desde ~/.conda/environments.txt
    if [ -f "${HOME}/.conda/environments.txt" ]; then
        while IFS= read -r env_dir; do
            if [ -x "${env_dir}/bin/python" ]; then
                CANDIDATES+=("${env_dir}/bin/python")
            fi
        done < "${HOME}/.conda/environments.txt"
    fi

    # Rutas comunes y entornos virtuales
    shopt -s nullglob
    CANDIDATES+=(
        "${REPO_DIR}/.venv/bin/python"
        "${REPO_DIR}/venv/bin/python"
        "${HOME}/miniconda3/envs/*/bin/python"
        "${HOME}/anaconda3/envs/*/bin/python"
        "${HOME}/.conda/envs/*/bin/python"
        "/opt/conda/envs/*/bin/python"
        /mnt/NAS/Tomas/results/*env*/bin/python
        /mnt/NAS/Tomas/*env*/bin/python3
        /mnt/NAS/Tomas/.venv/bin/python
        "${HOME}/.venv/bin/python"
        "${HOME}/miniconda3/bin/python"
        "${HOME}/anaconda3/bin/python"
    )
    shopt -u nullglob

    CANDIDATES+=(
        "$(command -v python3 2>/dev/null || true)"
        "$(command -v python 2>/dev/null || true)"
    )

    # Probar candidatos: se busca uno que tenga matplotlib y scipy
    for CAND in "${CANDIDATES[@]}"; do
        if [ -n "${CAND}" ] && [ -x "${CAND}" ]; then
            if "${CAND}" -c "import matplotlib, scipy" >/dev/null 2>&1; then
                PYTHON="${CAND}"
                break
            fi
        fi
    done
fi

# 3. Validar si se encontró un intérprete adecuado
if [ -z "${PYTHON}" ]; then
    echo "=========================================================="
    echo "[ERROR] No se encontro ningun entorno de Python con 'matplotlib' y 'scipy'."
    echo ""
    echo "Para ver tus entornos de Conda en este servidor, corre:"
    echo "    conda env list"
    echo ""
    echo "Luego podes activar tu entorno y correr:"
    echo "    conda activate <nombre_del_entorno>"
    echo "    bash run_caliali_batch.sh"
    echo ""
    echo "O pasar el entorno directamente:"
    echo "    bash run_caliali_batch.sh <nombre_del_entorno>"
    echo "=========================================================="
    exit 1
fi

# Asegurar que h5py esté presente en el Python seleccionado
if ! "${PYTHON}" -c "import h5py" >/dev/null 2>&1; then
    echo ">> Instalando h5py en el entorno (${PYTHON})..."
    "${PYTHON}" -m pip install h5py || true
fi

echo "=========================================================="
echo "Usando intérprete: ${PYTHON}"
${PYTHON} -c "import sys; print('Python: ' + sys.version.split()[0] + ' (' + sys.executable + ')')"
${PYTHON} -c "import matplotlib; print('Matplotlib: ' + matplotlib.__version__)"
${PYTHON} -c "import scipy; print('scipy: ' + scipy.__version__)"
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
