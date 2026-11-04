# Miniscope-Data

Procesamiento y generación de mapas de actividad espacial en CA1.

---

## Estructura del repositorio

```text
├── core/
│   ├── caliali_loader.py                # Carga de trazas CaliAli (S, C_raw, C) desde HDF5/MATLAB
│   ├── data_loader.py                   # Carga de tracking (x, y, t) y metadatos de DataBase
│   └── place_fields.py                  # Estimador gaussiano continuo de rate maps
├── run_batch.py                         # Script principal de procesamiento y generación de figuras
├── run_caliali_batch.sh                 # Wrapper para ejecutar en el servidor/clúster
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

Los scripts leen la trayectoria desde `DataBase` y la actividad de calcio desde `GDi` (CaliAli):
- **Trayectoria:** `/mnt/NAS/Miniscopes/Reg_CA1/DataBase`
- **CaliAli:** `/mnt/NAS/Miniscopes/Reg_CA1/GDi`

### 1. Ejecución mediante el script en el cluster

El script `./run_caliali_batch.sh` detecta automáticamente el entorno de Python y organiza las salidas en `results/`:

```bash
# Correr todo el lote completo con señal S (actividad deconvolucionada)
./run_caliali_batch.sh

# Correr todo el lote con fluorescencia cruda continua (C_raw)
./run_caliali_batch.sh --signal C_raw

# Correr solo un animal y una sesión específica (ej. para probar)
./run_caliali_batch.sh --signal C_raw --animal R006 --session 2026_06_16_merged.mat

# Asignar un nombre/sufijo personalizado a la carpeta de salida
./run_caliali_batch.sh --signal S --run-name nueva_version

# Forzar la regeneración y sobrescribir mapas si el .mat cambió
./run_caliali_batch.sh --signal S --animal R005 --session 2026_06_18_merged.mat --overwrite
```

### 2. Ejecución directa con python

```bash
python run_batch.py \
    --data-dir /mnt/NAS/Miniscopes/Reg_CA1/DataBase \
    --caliali-dir /mnt/NAS/Miniscopes/Reg_CA1/GDi \
    --signal S \
    --animal R005 \
    --session 2026_06_18_merged.mat
```
