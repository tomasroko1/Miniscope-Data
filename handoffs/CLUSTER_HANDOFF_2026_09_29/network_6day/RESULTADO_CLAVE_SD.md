# Qué encontramos en CA1, en dos minutos

La extensión a las dos tareas y diez jornadas está en
[el resultado para el congreso](RESULTADO_CONGRESO_CA1.md).

La prueba más amplia ahora incluye cuatro jornadas XsS y controla posición,
velocidad y actividad global. En las diez jornadas, la huella de coactividad
SAMPLE–TEST sigue superando los dos nulos tras ese control; ver la
[figura 18](18_CA1_maps_vs_networks.png). El detalle de abajo registra primero
el análisis SD previo, de seis jornadas.

**Resultado defendible:** la huella de coactividad de CA1 se conserva de SAMPLE
a TEST en las seis jornadas de la tarea S–D: tres con VEH y tres con CNO. No
observamos el colapso general de esa organización que se esperaría si los
conjuntos neuronales se desarmaran por completo bajo CNO. **No encontramos una
señal específica de S bajo CNO que explique el déficit conductual.**

La [figura principal](17_SD_pairwise_network.png) muestra las seis jornadas y,
separado a la derecha, el contraste CNO–VEH de cada animal. Si tenés que elegir
una sola figura neuronal para presentar, usaría esta.

## Qué significa la huella

En cada fase calculamos qué pares de neuronas tienden a activarse juntos. La
pregunta es si la lista de pares más relacionados en SAMPLE sigue siendo
similar en TEST, usando las **mismas identidades mapeadas dentro del día**.
Una similitud positiva podría aparecer por azar o por fluctuaciones globales;
por eso repetimos el cálculo permutando identidades y desplazando en el tiempo
la actividad de cada célula. También retiramos el modo global de la población.

| Animal | VEH: similitud observada / corte nulo 95% | CNO: similitud observada / corte nulo 95% |
|---|---:|---:|
| R004 | 0,193 / 0,020 | 0,106 / 0,009 |
| R005 | 0,110 / 0,017 | 0,286 / 0,021 |
| R006 | 0,570 / 0,034 | 0,636 / 0,036 |

El corte de esta tabla es el de desplazamientos temporales; los nulos por
identidad también se superan. Con 199 simulaciones por jornada, el mínimo `p`
empírico unilateral es 0,005 en cada una, y `p=0,03` tras Holm para las seis
jornadas. La dirección se mantiene con ventanas temporales de 0,5, 1 y 2 s
y contando solo el **inicio** de cada racha de actividad. Los seis días son
observaciones de **tres animales**, no seis réplicas biológicas independientes.

**Interpretación:** hay una arquitectura de coactividad reconocible que
atraviesa SAMPLE→TEST incluso en los días etiquetados CNO. La conducta S–D
alterada bajo CNO, si corresponde a estos mismos animales, no se puede
describir simplemente como “CA1 perdió toda su organización”. La prueba no
establece que el animal identifique correctamente S, que cada ensemble
discreto persista, ni cuál es el mecanismo del déficit. El contraste de
tratamiento de la huella, referido a OF1–OF2, es heterogéneo: CNO–VEH vale
aproximadamente −0,08, +0,16 y −0,32 para R004, R005 y R006.

Hay un contraste útil a otra escala: el **par de mapas espaciales menos
parecido** entre las cuatro fases tiene un margen de estabilidad inferior al
HabL del mismo animal en las 10/10 jornadas con objetos disponibles. Ocurre
con mapas de eventos y con mapas continuos; el cambio promedio por animal en
el margen de eventos es −0,054, −0,041 y −0,047 (R004/R005/R006). Es
compatible con que la representación espacial cambie al atravesar fases con
objetos **mientras persiste una estructura de coactividad**. No atribuye ese
cambio a CNO: HabL y objetos son fechas/contextos distintos y no hay HabL
pareado bajo ambos tratamientos. La
[tabla original de márgenes por jornada](../place_cell_stability_four_phases/sensitivity_cross200/strict_fourphase_six_pair_null_margin_by_day.csv)
permite revisar los 10 casos.

## Lo nuevo que probamos después

Medimos episodios de coactivación por encima de un corte nulo en regiones neuronales fijas
alrededor de los sitios S y D. En cada fase, el corte fue el percentil 99 de
100 nulos de desplazamientos circulares independientes por célula. Se contó
un máximo por episodio, separado al menos 0,5 s del siguiente. La tasa de
picos se dividió por segundos válidos **solo para normalizar la señal neural**;
no calculamos preferencia, permanencia ni conducta de acercamiento a objetos.
El índice fue el cambio TEST–SAMPLE en S menos D, restando la diferencia
espacial entre sitios virtuales de OF1.

| Animal | CNO−VEH, índice S−D de picos/s |
|---|---:|
| R004 | −0,140 |
| R005 | +0,288 |
| R006 | +0,020 |

La [figura de picos locales](SD_local_coactivity/paired_S_D_coactivity.png)
deja claro que los signos son mixtos. Cambiar la escala de cámara a 9 u 11
px/cm o ampliar la corona de 3–10 a 3–12 cm altera algunas direcciones. No
es un efecto S–D de CNO presentable. Los [datos por jornada](SD_local_coactivity/day_site_metrics.csv),
[contrastes pareados](SD_local_coactivity/paired_treatment_contrasts.csv) y
[protocolo](SD_local_coactivity/protocol.json) permiten revisarlo.

Vimos además una posible mayor sincronía global durante SAMPLE/TEST bajo CNO
en una definición de eventos. Igualamos el número de células entre fechas y
controlamos por bandas de velocidad. **No sobrevivió de modo uniforme** al
contar solo onsets o al binarizar la derivada de calcio con el criterio
Carrillo-Reid. Por haberla detectado después de mirar varias medidas, queda
como exploración, no como afirmación de “CNO aumenta la sincronía”. Los
[controles globales](SD_global_synchrony/paired_phase_contrasts.csv) y sus
sensibilidades [onsets](SD_global_synchrony/onsets/paired_phase_contrasts.csv)
y [derivada de calcio](SD_global_synchrony/calcium_derivative/paired_phase_contrasts.csv)
quedan guardados para auditoría.

## Frase lista para presentar

> En tres animales con jornadas S–D pareadas, la estructura de coactividad de
> CA1 entre SAMPLE y TEST fue mayor que dos controles aleatorios en las seis
> jornadas, tanto con vehículo como con CNO. No detectamos una reducción S–D
> específica bajo CNO que resistiera las definiciones y controles probados.
> El déficit conductual no parece corresponder a un colapso global evidente
> de esta huella de CA1, pero todavía no tenemos el índice conductual por
> jornada para vincular actividad neuronal y desempeño de cada animal.

## Límite que no debemos esconder

Las identidades celulares no se siguen entre días CNO y VEH, y tampoco está
confirmado que la gráfica conductual grupal incluya exactamente R004–R006.
La metadata no sitúa OF1 respecto de la administración de CNO. Por tanto,
podemos describir un correlato neuronal **compatible con que la organización
global de CA1 se conserva**, pero no una correlación individual conducta–CA1
ni causalidad. Una ausencia de reducción uniforme con `n=3` tampoco prueba
equivalencia entre tratamientos.

## Cómo reproducir esta comprobación local

Desde la raíz del proyecto: `python analyze_sd_local_coactivity.py` genera
los picos por sitio y sus controles; `python audit_sd_global_synchrony.py`
audita la sincronía global con el número de células igualado. Las dos
sensibilidades globales se repiten con `--representation onsets_S_gt_3SD` y
`--representation calcium_derivative_3SD`, usando sus respectivos
`--output-dir` para no sobrescribir los CSV principales. Las pruebas
sintéticas están en `test_sd_local_coactivity.py` y
`test_sd_global_synchrony.py`.
