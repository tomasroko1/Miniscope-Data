# Ensambles cross-day en el cluster

Este flujo usa el `mappings.csv` producido por `run_cross_registration_by_animal.py` y los `.mat` diarios de `DataBase/<ANIMAL>`. Empareja las adquisiciones de cada `.mat` con el registro por su ruta dia/adquisicion, y usa solo los `global_cell_id` presentes en las dos adquisiciones comparadas.

Por defecto corre cuatro comparaciones entre dias para fases equivalentes (OF1 a OF1, SAMPLE a SAMPLE, TEST a TEST, OF2 a OF2) y una comparacion TEST a OF2 entre dias. Cada ajuste concatena unicamente los dos bloques de esa comparacion. Stoixeion construye la matriz de similitud, asigna vectores a ensambles y guarda la reconstruccion temporal con el limite entre los dos bloques rotulado. No corre el ajuste de cuatro fases de una fecha ni enlaza ensambles por su numero E; la identidad comparable entre dias es el conjunto de `global_cell_id` del core.

## Comando

Desde la carpeta que contiene `coactivation_analysis/` y `core/`, corre en una sola linea. Reemplaza `cross_registration_TIMESTAMP` por la carpeta creada por cross-registration:

```bash
python3 coactivation_analysis/run_cross_day_ensembles.py --animal R005 --data-dir /mnt/NAS/Miniscopes/Reg_CA1/DataBase --mappings-csv /mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26/cross_registration_TIMESTAMP/mappings.csv --output-dir /mnt/NAS/Tomas/results/stoixeion/cross_day
```

La primera corrida usa 20 shuffles para los cortes de vectores y similitud. Se puede subir el minimo de celulas compartidas con `--minimum-cells 30`. El script informa que `.mat` omite si no tiene `act.mapping` verificable o si una adquisicion no aparece en `mappings.csv`.

## Salidas

- `cross_day_pairs.csv`: numero de celulas compartidas, estado, ensambles y directorio por comparacion.
- `cross_day_run.json`: parametros, pares incluidos y archivos omitidos.
- Una carpeta por par dia/fase: `carrillo_reid_figure.png`, `carrillo_reid_figure_vector_order.png`, matrices NPZ, resumen de ensambles y neuronas core.

El runner requiere archivos `.mat` mergeados de cuatro fases con `act.mapping` y `sess.sess_paths`. Los archivos `R005/2026_06_16_merged.mat` y `R005/2026_06_17_merged.mat` no tienen `act.mapping`; se omiten porque sus columnas locales no permiten traducir con seguridad a los `unit_id` del mapping cross-day.
