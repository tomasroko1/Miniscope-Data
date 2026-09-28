# Registro exploratorio de R004: HabL y SD

Esta corrida usa exactamente las 12 adquisiciones que figuran en `sess.sess_paths` de los `.mat` de R004 del 18, 19 y 21 de junio de 2026: HabL, SD CNO y SD VEH. El archivo `registration_manifests/R004_HabL_SD_CNO_VEH.txt` fija esos paths. No incluye HabC ni XsS.

## Qué porcentaje comparar

El mapping dentro de cada `.mat` solo vincula las cuatro fases de su fecha. Como referencia, las células presentes en las cuatro fases son:

| Fecha | Sesión | Células en 4/4 | Porcentaje según fase |
| --- | --- | ---: | ---: |
| 18/6 | HabL | 102 | 24,5–28,7 % |
| 19/6 | SD CNO | 169 | 29,3–34,6 % |
| 21/6 | SD VEH | 99 | 22,2–27,3 % |

Estos porcentajes intradía no son directamente comparables con exigir que una célula aparezca en las 12 adquisiciones. El nuevo reporte muestra porcentajes para cada par de adquisiciones, para cada par de fechas y el número de IDs presentes en las tres fechas. Los porcentajes por adquisici?n usan como denominador los centroides utilizables dentro del campo visual com?n a las 12 sesiones. Reducir la lista de sesiones puede ampliar el campo visual común y aumentar el número de células presentes en todas las sesiones elegidas; no garantiza que mejore la calidad o el porcentaje de cada par.

## Ejecución en el cluster

Activa un Python que pueda importar `minian`, `pandas`, `xarray`, `dask`, `skimage` y `matplotlib`. En la revision anterior, el Python del sistema y el entorno `miniscope-habl` no tenian `minian`. Si sigue faltando, usa el metodo de [instalacion con conda-forge recomendado por MiniAn](https://minian.readthedocs.io/en/stable/start_guide/install.html):

```bash
source /mnt/NAS/Tomas/miniconda311/bin/activate
conda create -y -n minian-registration
conda activate minian-registration
conda install -y -c conda-forge minian matplotlib
python -c "import minian, pandas, xarray, dask, skimage, matplotlib; print('ready')"
```

La opcion `--require-qc-figures` evita una corrida sin las figuras solicitadas. Desde cualquier directorio, corré en una sola línea:

```bash
python /mnt/NAS/Tomas/Miniscope-Data/run_cross_registration_by_animal.py --animal-root /mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26 --session-list /mnt/NAS/Tomas/Miniscope-Data/registration_manifests/R004_HabL_SD_CNO_VEH.txt --reference-mat-dir /mnt/NAS/Miniscopes/Reg_CA1/DataBase/R004 --require-qc-figures
```

La salida se crea en una carpeta nueva `cross_registration_YYYYMMDD_HHMMSS` debajo de `R004_M4_26`; el log imprime el path exacto. El script comprueba que estén las 12 carpetas del manifiesto antes de abrir los datos.

## Salidas principales

- `mappings.csv`, `mappings_matlab.csv` y `cross_registration_log.log`.
- `pairwise_match_summary.csv`: matches, porcentajes por adquisición, distancias entre centroides y área de superposición de cada par.
- `within_day_mapping_comparison.csv` y `reference_unit_id_coverage.csv`: compara los pares de unit_id de la nueva corrida con los mappings intradía de los tres `.mat`, con porcentajes y coincidencias exactas.
- `day_match_summary.csv` y `registration_scope.json`: coincidencias entre HabL, SD CNO y SD VEH; campo visual común a las 12 adquisiciones.
- `centroid_fov_overview.png`: proyecciones alineadas, centroides utilizables, centroides enlazados entre fechas y borde del campo visual común.
- `centroid_matches_*.png`: seis láminas (pares de fechas y comparaciones dentro de cada fecha) con proyecciones superpuestas, centroides, líneas de matches y bordes de las áreas compartidas.
- `centroid_match_pairs.csv.gz`: coordenadas y distancias de los pares dibujados.

El nombre de HabL indica posibles franjas en dos adquisiciones. Revisá sus paneles de alineación y los centroides antes de interpretar los matches entre fechas.
