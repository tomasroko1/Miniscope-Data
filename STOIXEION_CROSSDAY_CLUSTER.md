# Stoixeion cross-day en MATLAB

El registro de células entre sesiones se genera con MiniAn en Python y produce `mappings.csv`. La detección de ensembles se corre en MATLAB con `run_stoixeion_crossday.m` y la versión de Stoixeion con `stoixeion_exports.patch` aplicada. El runner MATLAB analiza pares de fechas usando únicamente `global_cell_id` presentes en ambas adquisiciones. Compara las cuatro posiciones de subsesión y agrega TEST→OF2 cuando ambas fases existen. No ejecuta la concatenación global de las cuatro fases de una fecha.

## Registro celular

Desde el entorno Python que tiene MiniAn y desde el repo:

```bash
python3 /mnt/NAS/Tomas/Miniscope-Data/run_cross_registration_by_animal.py --animal-root /mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26
```

El script guarda una carpeta `cross_registration_...` debajo de la carpeta del animal. `matplotlib` es opcional: si no está instalado, el script registra un aviso y omite las figuras de control, pero sigue con el mapping. Esta etapa Python solo genera el mapping.

## Ensambles en MATLAB

Abrí MATLAB y pegá la siguiente llamada en la **ventana de comandos de MATLAB**, no en Bash. Usa la ruta real a `mappings.csv`:

```matlab
run_stoixeion_crossday('/mnt/NAS/Miniscopes/Reg_CA1/DataBase','/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/cross_registration_TIMESTAMP/mappings.csv','/mnt/NAS/Tomas/results/stoixeion/cross_day','R005','/mnt/NAS/Tomas/Stoixeion')
```

El runner guarda `cross_day_pairs.csv`, las figuras nativas de Stoixeion, `cross_day_ensemble_timeline.png` y `core_neurons.csv` por cada comparación. Las sesiones simples sin `act.mapping` se omiten porque no se pueden traducir con seguridad a las identidades del registro cross-day.
