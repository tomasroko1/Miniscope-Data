# Un resultado neuronal presentable de CA1

**Figura para insertar:** [mapas y red, diez jornadas](18_CA1_maps_vs_networks.png).
Cada fila es una jornada con objetos. Hay tres animales y diez jornadas con
mapping inequívoco: seis S–D y cuatro XsS, cinco VEH y cinco CNO.

## Lo que medimos

**Mapas individuales, panel izquierdo.** Para cada célula seguimos su mapa
entre OF1, SAMPLE, TEST y OF2 del mismo día. De las seis comparaciones entre
fases tomamos la menos estable, descontando el nivel esperado por
desplazamientos temporales aleatorios. Resumimos las células de cada jornada
y comparamos ese valor con el HabL del **mismo animal, registrado otro día**.
El margen es menor en **10/10 jornadas con objetos**, tanto con frames-evento
`S>3 SD` como con actividad `S` continua. Una jornada no equivale a un animal
independiente: las diez jornadas pertenecen a tres animales.

**Red de coactividad, panel derecho.** Para cada jornada usamos solo células
mapping identificadas en las cuatro fases. Contamos una vez el inicio de cada
racha `S>3 SD`, reunimos la actividad en ventanas de 1 s y calculamos qué
pares de células tienden a activarse juntos. La estructura de pares de SAMPLE
reaparece en TEST en **10/10 jornadas**. En todas supera tanto la permutación
de identidades como el desplazamiento temporal independiente de cada célula,
incluso después de predecir y quitar, para cada neurona y fase, la actividad
asociada con posición, velocidad y fluctuación global. El control de posición
explica una mediana de 6–12% de la variación por fase. Son 199 nulos de cada
clase por jornada: `p` empírico unilateral = 0,005 en cada una; corrección
Holm entre diez jornadas = 0,05. Al contar todos los frames-evento en vez de
solo onsets, también se reproduce en 10/10 tras el mismo control.

El cociente visual del panel derecho es similitud observada dividida por el
corte 95% **más alto** de los dos nulos, después del control de posición y
velocidad. Vale al menos 2,37. Es una ayuda para ver la separación del nulo,
**no** un tamaño de efecto farmacológico. El modelo de posición es una
corrección razonable, no una descripción exhaustiva de todos los movimientos.

## Qué decir en el congreso

> Durante las tareas con objetos, los mapas espaciales de células individuales
> muestran menor estabilidad entre fases que en HabL del mismo animal. A la
> vez, la organización de pares de neuronas que coactivan se conserva de
> SAMPLE a TEST, tanto con vehículo como con CNO. Esto sugiere que CA1
> reconfigura representaciones sin perder por completo su estructura de red.

Esa es una **observación neuronal** que va más allá de mostrar mapas bonitos.
La comparación con HabL sigue siendo descriptiva porque corresponde a otra
fecha. La persistencia de red se prueba dentro de cada jornada; no implica
que el animal recuerde bien el objeto desplazado ni que cada ensemble de
Carrillo-Reid permanezca idéntico. La hipótesis más interesante para discutir
es que el fallo conductual bajo CNO podría afectar **qué** información se
separa o se lee de una red todavía organizada, pero eso **no quedó demostrado**.

## Qué no afirmar

- No hay una diferencia CNO−VEH en CA1 que sea uniforme en los tres animales
  para la discriminación neuronal específica de S, la estabilidad de mapas,
  los cores Carrillo-Reid ni las medidas alternativas que probamos.
- La presentación conductual no da un índice S–D por jornada ni confirma que
  su grupo conductual incluya exactamente R004–R006. No calculamos una
  correlación individual conducta–CA1.
- El día XsS-VEH de R004 (30 de junio) tiene 540 columnas de actividad en OF1
  y solo 539 IDs mapping finitos; no asignamos una identidad al azar. Se
  excluye de las pruebas que requieren seguir la misma célula entre fases.
  Sus mapas se pueden evaluar fase por fase, pero no corregir esa asignación
  con la información disponible.
- No hicimos análisis conductual de cercanía ni tiempo junto a objetos.

## Archivos para auditoría

- [Diez redes: onsets y nulos](All_object_networks/all_object_networks_onsets.csv)
  y [sensibilidad con frames-evento](All_object_networks/all_object_networks_frames.csv).
- [Red tras controlar posición y velocidad: onsets](All_object_networks/position_controlled_networks_onsets.csv)
  y [frames](All_object_networks/position_controlled_networks_frames.csv).
- [Tabla que une mapas y redes por jornada](CA1_context_network_joint_evidence.csv).
- [Margen original de mapas](../place_cell_stability_four_phases/sensitivity_cross200/strict_fourphase_six_pair_null_margin_by_day.csv).
- Código: `audit_all_object_networks.py`, `audit_position_controlled_networks.py`
  y `make_ca1_context_network_figure.py` en la raíz del proyecto. El test
  sintético de identidad/permutación se ejecuta con
  `python audit_sd_pairwise_network.py --self-test`.

## Candidatas espaciales tipo place cell en la tarea SD

También comparé HabL con SAMPLE y TEST usando el criterio ya calibrado con los
nulos de HabL: una célula cuenta como **candidata espacial** si su información
espacial y su correlación split-half superan ambos cortes HabL (`q95`:
2,132 bits/evento y `r=0,267`). Es una clasificación por fase, no una
afirmación de que cada candidata sea una place cell confirmada.

En SAMPLE la fracción de candidatas `q95` fue menor con CNO en los tres
animales:

| Animal | VEH | CNO | CNO−VEH |
|---|---:|---:|---:|
| R004 | 5,05% | 3,55% | −1,50 pp |
| R005 | 11,03% | 8,48% | −2,55 pp |
| R006 | 4,62% | 1,20% | −3,42 pp |
| Media de animales | 6,90% | 4,41% | −2,49 pp |

Al subir el corte a `q99`, la dirección ya no es uniforme: R004 baja de 2,02%
a 0,59%, R005 sube de 2,21% a 2,42% y R006 queda en 0% en ambos días. En TEST
el corte `q95` es menor con CNO en R004 y R006, y prácticamente igual en R005
(7,07%→1,78%; 6,62%→6,67%; 4,62%→3,61%, VEH→CNO). La información espacial
mediana y la estabilidad media de mapas también varían por animal. Por eso
presentaría el resultado de SAMPLE como una **tendencia piloto en candidatas
espaciales**, no como pérdida confirmada de mapas por CNO.

La [figura de comparación tipo place cell](19_SD_place_like_comparison.png)
muestra tres filas, una por animal, y un punto por neurona. Compara HabL OF1
con VEH y CNO en SAMPLE y en TEST. El orden agrupa las dos condiciones de
SAMPLE y después las de TEST para facilitar la comparación por fase; no es
una secuencia entre días. La columna izquierda mide información espacial; la
derecha mide fiabilidad split-half. El
[resumen por animal](SD_place_like_by_animal.csv) conserva denominadores y
métricas, incluidos los porcentajes que pasan `q95` y `q99`. Los puntos
verdes pasan ambos cortes q95. Las líneas horizontales son cortes de
clasificación por neurona, no pruebas estadísticas CNO–VEH. Las cajas
describen neuronas anidadas en tres animales; las observaciones celulares no
cuentan como réplicas biológicas independientes. TEST queda en la tabla CSV,
no en esta figura.
