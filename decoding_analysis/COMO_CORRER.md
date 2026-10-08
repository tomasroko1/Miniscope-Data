# Decodificación y células de objetos

Ejecutar desde la raíz del repositorio. No cambia `core/` ni necesita dependencias
nuevas. La lógica pequeña de nombres de fases y categorías replica `run_batch.py`
para no importar su runner de figuras. Ni el loader ni el código revisado guardan
coordenadas de objetos: el análisis las obtiene exclusivamente del CSV.

## Entorno

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.txt
export MPLBACKEND=Agg
export MPLCONFIGDIR=/tmp/mpl
export OPENBLAS_NUM_THREADS=1
export MINISCOPE_DATA_DIR=/ruta/a/DataBase
```

La instalación es sólo para preparar un entorno nuevo con acceso a los paquetes;
en este workspace ya existe `.venv` y no se instala nada. El valor de datos por
defecto es `/mnt/NAS/Miniscopes/Reg_CA1/DataBase`. Estructura:
`DATA_DIR/animal/archivo.mat`. `--data-dir` tiene prioridad sobre la variable.
Limitar BLAS a un hilo evita sobrecostes en las numerosas multiplicaciones pequeñas.

## Comandos

```bash
.venv/bin/python -m decoding_analysis.run_analysis decode --data-dir "$MINISCOPE_DATA_DIR" --out-dir results --time-bin 0.5 --decoder-bin 2 --sigma 4 --k-folds 5 --fold-gap-s 2 --min-speed 2 --cells common --n-shuffles 100 --save-predictions
.venv/bin/python -m decoding_analysis.run_analysis objects --objects-csv decoding_analysis/object_coords.csv --radius 5 --min-speed 0 --min-events 10 --n-shuffles 500 --alpha 0.05 --n-examples 6
.venv/bin/python -m decoding_analysis.run_analysis all --objects-csv /ruta/coordenadas_reales.csv --seed 0
```

Filtros opcionales: `--animal Mouse1 --session T1_SD_CNO.mat`. `--session` es el
nombre de archivo con o sin `.mat`, no `day_name`. `--activity events|S` (default events),
`--threshold-sd 3`, `--prior uniform|occupancy`, `--equalize-duration`.
`--decode-shuffles` y `--object-shuffles` permiten valores distintos en `all`;
`--n-shuffles` sobrescribe ambos. En `all`, sin sobrescribir, se usan 100/500.
`--min-speed` y `--min-events` también sobrescriben ambos análisis: si no se
especifican, decoding usa 2 cm/s y 1 evento en training; objetos usa 0 y 10.
`--min-object-occupancy` fija ocupación mínima tanto dentro como fuera (2 s).
`--arena-tolerance-cm 2` acepta posiciones entre radio 25 y 27 cm, proyectándolas
radialmente a 25 cm menos un epsilon; `0` desactiva esta tolerancia. Por encima
se excluyen. Si más del 5% de todos los frames de una fase exceden el radio con
tolerancia, se imprime **¡¡¡ ARENA !!!**: revisar centrado en (0,0) y unidades cm.
Cada ejecución guarda todos los argumentos CLI y los defaults efectivos de filtros
y shuffles en `<out-dir>/run_config.json`; usar directorios separados para conservar
la procedencia de distintas ejecuciones.
Errores se registran en `analysis_failures.csv`, el runner continúa con las demás
sesiones/análisis y termina con código 1 si hubo fallos. `--verbose` imprime traceback.

## Coordenadas reales

Editar `decoding_analysis/object_coords.csv` o pasar otro CSV. Las filas incluidas
son **PLACEHOLDER**, no datos experimentales: dos pares a 15/25 cm (SD), 5/15 cm
(XsS), y dos objetos TEST ilustrativos. Se emite una advertencia explícita al
utilizarlas y se registra `coords_are_placeholder=True`. No interpretar esos
resultados como evidencia sobre objetos reales.

Columnas obligatorias:

| Columna | Significado |
|---|---|
| session_type | SD o XsS |
| day_name | etiqueta exacta `sess.day_name`, vacía = todas |
| animal | animal exacto, vacío = todos |
| phase | SAMPLE o TEST |
| object_id | identidad única dentro de fase, p.ej. A1..A4 / T1,T2 |
| x_cm, y_cm | centímetros, mismo marco que `pos.x/y`, centro arena (0,0) |
| match_id | en TEST: `object_id` SAMPLE correspondiente; vacío = sin correspondencia |
| is_placeholder | 1 para coordenadas ilustrativas, 0 para medidas reales |
| note | texto explicativo |

Cada fase elige **un conjunto completo**: animal+day_name > day_name > sólo animal
> session_type sin overrides. No combina coordenadas de distintos niveles.
Definir todas las filas de la fase en cada override. Se permiten comentarios `#`.
Si se mezclan fases a distintos niveles, cada una se selecciona independientemente.
Identidades duplicadas, correspondencias múltiples y posiciones fuera de arena son
errores. `match_id` controla exclusivamente stayed (distancia <2 cm), moved (>=2 cm)
o removed (ninguna correspondencia). No se presupone qué objetos retiene TEST.
OF1/OF2 utilizan coordenadas SAMPLE como pseudo-objetos.

## Métodos y supuestos

* Arena circular de radio 25 cm, centrada en (0,0). Posición y velocidad en cm y
  cm/s; `pos.v` se utiliza directamente. Si falta se calcula derivada de x/y/t,
  suavizada a 0.25 s. Para calcular esa velocidad se interpolan posiciones faltantes;
  los frames originales con x/y NaN siguen inválidos. Intervalos mayores que
  max(0.1 s, 3 veces la mediana de dt) se consideran huecos: el frame anterior se
  invalida y su dt se limita a la mediana para no atribuir el hueco a ocupación.
  Si falta `act.t`, el loader usa 0.05 s (20 Hz).
* Decoding: bins de 0.5 s; actividad sumada, posición y velocidad medias.
  Excluye bins vacíos, incompletos/con frames inválidos sin descartar la fase; sólo velocidad >=2 cm/s por
  defecto en entrenamiento y evaluación. k=5 bloques contiguos del tiempo
  original, con gap 2 s por lado. No concatena periodos de movimiento.
* Poisson naive Bayes en log espacio, MAP y error de media posterior. El modo S
  trata amplitudes positivas como pseudocuentas: la hipótesis Poisson es aproximada.
  Mapas: histogramas de ocupación/actividad de **training**, suavizados con sigma
  4 cm, normalización por ocupación y piso 1e-6. Es la versión discretizada del
  estimador gaussiano continuo del laboratorio. Default espacial 2 cm para ahorrar
  cálculo en CV/shuffles; `--decoder-bin 1` reproduce la resolución de laboratorio.
  El grid cubre exactamente [-25,25]; si 50 no es múltiplo del bin solicitado,
  el espaciado efectivo es 50/ceil(50/bin). Prior uniforme por defecto; occupancy
  usa ocupación suavizada de training. Azar: una posición aleatoria por bin de test,
  uniforme o según histograma de ocupación cruda de training según el prior.
* events reproduce `threshold_event_frames` de forma vectorizada por columnas:
  S finita > 3*SD de los valores fuente finitos (ddof=1); NaN/Inf produce 0. En decoding, SD y
  selección por actividad se calculan sólo en training de cada fold (también en
  cada shuffle); en objetos se calcula por fase. No normaliza con bins de test.
  `n_cells` indica población disponible; `n_cells_used` es la mediana entre folds
  de células que pasan min-events (cada fold pesa igual); `n_cells_min_fold`/predictions indican
  población que supera `--min-events` en cada fold.
* common usa `get_common_neurons` y `mapping_columns`. all usa cada población local,
  conservando global IDs si hay mapping válido. Mapping ambiguo/ausente se anuncia y
  cae a columnas locales por fase; entonces no compara identidades entre fases ni
  calcula desplazamiento o flags OF cruzados. IDs son zero-based; `local_column`
  siempre refiere la columna original de S. Sesiones simples tienen una sola fase.
* Equalize-duration conserva el mismo número de bins de movimiento en todas las
  fases (primeros N elegibles). CV/shifts mantienen el timeline completo, aunque
  sólo esos bins participan en mapas y métricas.
* Shuffles: shifts circulares comunes a la población, mínimo 30 s en ambos sentidos,
  sobre timeline completo, no sobre bins concatenados después de filtrar velocidad.
  Decoding desplaza actividad bineada y vuelve a entrenar cada fold. Reutiliza
  resultados de shifts repetidos (conservando multiplicidad en el nulo), y en modo
  S reutiliza actividad bineada independiente del fold/shift. Objetos desplaza
  actividad por frames, con FFT vectorizada. Con shuffles requiere duración >60 s; se asume
  muestreo aproximadamente regular (shifts en índices, no interpolación temporal).
  Un shift común preserva correlaciones entre células. p decoding = (1+#nulo<=obs)/(N+1),
  p objetos = (1+#nulo>=obs)/(N+1). z decoding = (obs-media_nulo)/SD_nulo: negativo
  significa error menor. Con N=0 las estadísticas del nulo son NaN.
* Objetos: zona unión de círculos (sin contar dos veces solapamientos), tasa = suma
  actividad/ocupación. Índice (in-out)/(in+out). Cada objeto tiene además tasa propia,
  que puede compartir frames con otro objeto próximo. Preferido = mayor tasa.
  Criterio solicitado: p<alpha e índice>0; se publica q BH por fase y también
  `is_object_cell_fdr` (q<alpha). Células no elegibles: p/q NaN y motivo explícito.
  `sig_in_OF_control` señala enriquecimiento significativo en cualquier OF, para
  distinguir posibles campos de lugar coincidentes. `is_object_cell_strict` requiere
  `is_object_cell=True` y `sig_in_OF_control=False`: excluye significancia en
  cualquiera de los OF. Si falta mapping, una fase OF o la célula no tiene ambos
  controles y ninguno es significativo, el flag OF es NaN y el estricto False.
  Un control presente con actividad insuficiente no cuenta como significativo.
  **El enriquecimiento de
  proximidad por sí solo no demuestra selectividad causal al objeto.**
* Desplazamiento usa células significativas en SAMPLE con mapping alineable;
  si falta la célula en TEST, conserva clasificación y deja tasas NaN.
  Reporta tasas TEST en lugares viejo/nuevo, follow_index sólo para moved, correlación
  mapas SAMPLE/TEST, distancia del pico TEST a ambas ubicaciones y significancia TEST.
  Figuras/mapas de desplazamiento usan `gaussian_rate_maps_matlab` a 1 cm/sigma 4,
  respetando el filtro de velocidad. No se usan estos mapas para entrenar decoding.
* Comparación objetos−OF: media de errores medianos SAMPLE/TEST menos OF1/OF2 por
  sesión. Se añaden `objects_shuffle_z`, `of_shuffle_z` (medias respectivas) y
  `delta_shuffle_z` = objetos menos OF. Un delta z negativo favorece objetos
  respecto a su propio nulo; complementa el error crudo dependiente de ocupación
  y duración, sin garantizar comparabilidad completa. Wilcoxon bilateral si n>=3, global y por tratamiento; de otro modo NaN con
  nota. Sesiones repetidas de un animal no son réplicas independientes: tabla
  descriptiva solicitada, no modelo jerárquico de inferencia por animal.

## Salidas

`results/decoding/<animal>/<session>/decoding_<session>.csv` y
`results/decoding/decoding_summary.csv`: metadata animal/session_file/day_name/task/
treatment/phase/has_objects, población, bins y tiempo, filtros, error mediano/medio,
error de media posterior, media/mediana/p5 del nulo, p/z y azar.
`frac_frames_projected` y `frac_frames_outside_arena` se calculan por fase sobre
**todos los frames**, incluidos los que se invalidan por otras causas; también
aparecen en las tablas de células de objetos. Fuera de arena significa radio
mayor que 25+tolerancia; x/y no finitos se invalidan, no se cuentan como fuera.
`shuffle_<phase>.csv` guarda cada mediana nula; `predictions_<phase>.csv` opcional
guarda bin, tiempo, fold, posiciones real/MAP, velocidad, células y errores.
`decoding_objects_vs_of.csv` y `decoding_objects_vs_of_wilcoxon.csv` guardan deltas y
pruebas. `figures/` contiene curvas por sesión con nulo/azar, histogramas por sesión,
deltas por tratamiento, y HabL separado.

`results/object_cells/<animal>/<session>/object_cells_<session>.csv` y resumen:
metadata, global_id/id_scope/local_column, n_events (suma de actividad en modo S),
tasas in/out/per objeto, índice, preferido, p/q, flags p/FDR/OF/estricto, placeholder, radio,
ocupación, motivo de exclusión y filtro. Displacement por sesión y
`object_displacement_summary.csv` contienen las métricas anteriores.
`object_cell_counts.csv` por sesión y combinado contiene población/elegibles,
conteo p/FDR/estricto (`n_object_cells_strict`) y fracciones
(`fraction_object_cells`, `fraction_object_cells_strict`) (denominador: todas las células presentes) por fase/tratamiento.
`example_<ID>.png` muestra trayectoria/actividad, mapa y nulo por fase para top N
SAMPLE (estrictas primero, después crudas; luego p ascendente e índice descendente). No enlaza columnas locales entre fases.

## Prueba sintética

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/mpl OPENBLAS_NUM_THREADS=1 .venv/bin/python -m decoding_analysis.make_synthetic --out-dir /tmp/synth_data --seed 0
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/mpl OPENBLAS_NUM_THREADS=1 .venv/bin/python -m decoding_analysis.run_analysis all --data-dir /tmp/synth_data --out-dir /tmp/synth_results --objects-csv /tmp/synth_data/object_coords_synthetic.csv
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/mpl OPENBLAS_NUM_THREADS=1 .venv/bin/python -m decoding_analysis.check_synthetic --data-dir /tmp/synth_data --results-dir /tmp/synth_results
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/mpl OPENBLAS_NUM_THREADS=1 .venv/bin/python -m unittest decoding_analysis.test_analysis
```

Defaults: 2 animales, 5 sesiones por animal, 240 s/fase, 20 Hz, 80 células.
Opciones del generador: `--duration-s`, `--n-cells`, `--n-animals`, `--seed`.
Genera MATLAB cells con arrays object 1-D incluso para fases de igual longitud,
IDs locales permutados y ausencias NaN. Comprueba tipo merged y alineación con
el loader real. Produce ground_truth.csv y CSV de coordenadas sin placeholders.
Ocho células de objetos por sesión task: cuatro stayed y cuatro moved; las demás
son campos de lugar estables (~60%) o ruido. TEST conserva A1 y mueve A3.
Para iterar: `--n-shuffles 20 --n-examples 1` en otro out-dir; 20 shifts sólo dan
p mínimo 1/21, por lo que q FDR puede no detectar células plantadas.
El checker exige todos los archivos/fases, errores <85% del p5 shuffle y azar,
recall mínimo 70%, falsos positivos crudos <=20% (campos de lugar coincidentes incluidos), falsos positivos estrictos <=10%, imprime recall estricto, y exige clasificación
correcta de moved/stayed, con follow_index positivo en >=80% de moved.
