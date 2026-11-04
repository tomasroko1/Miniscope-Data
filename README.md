# Miniscope-Data

Pipeline para el procesamiento y generación de mapas de actividad espacial (*rate maps* / *place cells*) en hipocampo (CA1), integrando trazas unificadas de **CaliAli** y datos de comportamiento/tracking animal de **DataBase**.

---

## Estructura del repositorio

```text
├── core/
│   ├── caliali_loader.py                # Carga de trazas CaliAli (S, C_raw, C) desde HDF5/MATLAB
│   ├── data_loader.py                   # Carga de tracking (x, y, t) y metadatos de DataBase
│   └── place_fields.py                  # Estimador gaussiano continuo de rate maps
├── run_batch.py                         # Script principal de procesamiento y generación de figuras
├── run_caliali_batch.sh                 # Wrapper para ejecutar en el servidor/clúster
├── generar_diagnostico_dlc_excluidos.py # Control de calidad para sesiones con saltos de DeepLabCut
└── requirements.txt                     # Dependencias mínimas de Python
```

---

## Instalación

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

---

## Uso

Los scripts leen el tracking desde `DataBase` y las trazas unificadas desde `GDi` (CaliAli):
- **Tracking:** `/mnt/NAS/Miniscopes/Reg_CA1/DataBase`
- **CaliAli:** `/mnt/NAS/Miniscopes/Reg_CA1/GDi`

### 1. Ejecución mediante el script bash (Recomendado en clúster)

El script `./run_caliali_batch.sh` detecta automáticamente el entorno de Python y organiza las salidas en `results/`:

```bash
# Correr todo el lote completo con señal S (espigas deconvolucionadas)
./run_caliali_batch.sh

# Correr todo el lote con fluorescencia cruda continua (C_raw rectificada en 0)
./run_caliali_batch.sh --signal C_raw

# Correr solo un animal y una sesión específica (ej. para probar)
./run_caliali_batch.sh --signal C_raw --animal R006 --session 2026_06_16_merged.mat

# Asignar un nombre/sufijo personalizado a la carpeta de salida
./run_caliali_batch.sh --signal S --run-name nueva_version

# Forzar la regeneración y sobrescribir mapas si el .mat cambió
./run_caliali_batch.sh --signal S --animal R005 --session 2026_06_18_merged.mat --overwrite
```

### 2. Ejecución directa con Python

```bash
python run_batch.py \
    --data-dir /mnt/NAS/Miniscopes/Reg_CA1/DataBase \
    --caliali-dir /mnt/NAS/Miniscopes/Reg_CA1/GDi \
    --signal S \
    --animal R005 \
    --session 2026_06_18_merged.mat
```

---

## Parámetros clave

| Parámetro | Opciones | Descripción |
| :--- | :--- | :--- |
| `--signal` | `S`, `C_raw`, `C` | Tipo de señal: `S` (espigas deconvolucionadas, default), `C_raw` (fluorescencia cruda rectificada en 0) o `C` (traza continua denoised). |
| `--animal` | `R004`, `R005`, `R006` | Filtra por un animal específico. |
| `--session` | Archivo `.mat` | Filtra por un archivo específico (ej: `2026_06_19_merged.mat`). |
| `--run-name` | Texto | Crea una carpeta con nombre personalizado: `results/Mapas_por_tipo_caliali_<signal>_<run-name>`. |
| `--overwrite` | Flag | Regenera las figuras existentes en vez de omitirlas. |

---

## Resultados y exclusiones

- Las figuras se guardan jerárquicamente por animal y condición (ej. `results/Mapas_por_tipo_caliali/R005/Objetos/SD/VEH/2026_06_19_merged/neuron_0000.png`).
- Cada figura alinea las fases de la sesión (2 filas × N fases: trayectoria arriba y mapa de actividad abajo), manteniendo la misma escala visual para poder comparar fases.
- **Sesiones excluidas por DeepLabCut:** `R004-XsS-CNO` y `R005-XsS-VEH` poseen saltos espurios de tracking. Son omitidas automáticamente del cálculo de mapas y documentadas con `trayectoria_anomala_deeplabcut.png` y `LEEME_EXCLUSION_DLC.txt`.
