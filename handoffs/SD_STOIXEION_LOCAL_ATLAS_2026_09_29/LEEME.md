# Atlas Stoixeion S–D reconstruido localmente

Abrí `index.html`. La lámina `comparacion_seis_jornadas.png` pone VEH y CNO lado a lado para cada animal, con una escala de color compartida **dentro de ese animal**. Además hay una figura alta con **todos** los factores globales de cada día, un panel tipo Carrillo-Reid, las listas de células miembro y core, los frames asignados y la matriz binaria de coseno muestreada. `factor_quality.csv` añade banderas **descriptivas** para cores de menos de tres células, membresía mayor al 80 % de las células mapeadas y asignaciones en menos del 1 % de los frames; no descarta factores ni constituye una prueba estadística.

## Procedencia y cálculo

- **Identidad de los factores y frames asignados:** ajuste MATLAB Stoixeion ya exportado en `results/stoixeion_fast_clean/`; cuatro fases concatenadas por jornada. Este script **no volvió a ajustar la SVD global**: para R004 CNO la matriz de 32.877 vectores ocuparía más de 8 GB solo como arreglo `float64`, antes de las copias y SVD, mientras que esta PC tiene 8 GB de RAM y no tiene MATLAB.
- **Cálculo local nuevo:** lectura de los seis `.mat` de `data/`, mapping de las células presentes en las cuatro fases, umbrales de `S` exportados (`3×SD(S)` global por célula), reconstrucción de los vectores de alta actividad, validación de cada frame asignado y de los conteos por factor/fase, curvas en ventanas de 10 s, métricas de células core y matriz de coseno TF-IDF binarizada de una muestra temporal equilibrada entre fases. `QC_ALL_DAYS.json` comprueba las seis jornadas y registra SHA-256 de cada `.mat`.
- **Matriz M:** hasta 110 vectores de alta actividad por fase; TF-IDF usa todos los vectores significativos del día, `scut` el umbral del ajuste exportado. La matriz se muestra **antes** de los dos filtros Jaccard y la SVD. No es la matriz completa.
- **Eje temporal:** `OF1 → SAMPLE → TEST → OF2`; cada raya es un frame asignado; la curva muestra porcentaje de **todos** los frames asignados en cada ventana de 10 s. La escala vertical cambia por factor y está impresa en cada fila. `phase_metrics.csv` permite comparar valores exactos entre fases.
- **IDs de células:** `mapped_cell_id` del `.mat`, base cero, válido entre fases del mismo día. No identifica la misma célula entre días VEH/CNO.

## Lectura científica

Esto permite ver cuándo se expresa cada factor a lo largo de una jornada y descargar la lista de sus células. Un factor global fue definido usando **las cuatro fases**; su presencia en SAMPLE y TEST no es una prueba independiente de reactivación. Cada marca es un frame, no un episodio independiente. Los paneles describen seis jornadas anidadas en tres animales y no demuestran por sí solos un efecto de CNO, codificación del objeto desplazado ni una explicación del comportamiento.

## Reproducción

Desde la raíz de `revision-datos`:

```powershell
python coactivation_analysis/build_sd_stoixeion_atlas_local.py --data-dir data --exports-dir results/stoixeion_fast_clean --out results/SD_stoixeion_atlas_local_20260929
```
