# Qué decir mañana sobre Carrillo-Reid / Stoixeion

**Decisión para el póster:** estas figuras documentan candidatos de ensemble y continuidad de membresía celular, pero **no identifican un efecto de CNO ni explican la conducta**. El 5/6 es exploratorio y no debe ocupar el panel neuronal principal. Para ese panel, usar la huella de coactividad **por pares** de 6/6 jornadas S–D, rotulada como análisis separado.

## Respuesta de 25 segundos

> «Usamos Stoixeion para detectar patrones recurrentes de coactividad de CA1 sin indicar al algoritmo cuándo había objetos. Lo corrimos **por separado en cada fase** y seguimos entre fases las células *core* de los patrones detectados. En **5 de las 6 jornadas S–D** hubo al menos un par de cores SAMPLE–TEST con más células compartidas que las esperadas al permutar identidades, después de corregir las comparaciones múltiples. En R006 CNO, dos cores pueden seguirse durante OF1, SAMPLE, TEST y OF2. Esto muestra continuidad de la composición de esos grupos **dentro del día**; no demuestra todavía qué objeto codifican ni un efecto uniforme de CNO».

## Figuras de apoyo si preguntan por Stoixeion

1. **Resultado de las seis jornadas:** [solapamiento core SAMPLE–TEST](SD_core_overlap_6_days.png). Cinco jornadas tienen al menos un par con `q≤0,05`; la excepción es **R004 S–D CNO**. Son tres animales con un día VEH y uno CNO cada uno. El `q` procede de **999 permutaciones**, con Benjamini–Hochberg sobre los 176 pares SAMPLE–TEST posibles de diez jornadas con objetos. [Números de la figura](SD_core_overlap_6_days.csv).
2. **Ejemplo temporal:** [dos cores de R006 S–D CNO a lo largo de las cuatro fases](R006_CNO_two_core_timecourses.png). La curva es la fracción de frames asignados en ventanas de 10 s; las marcas altas son frames individuales. Los nombres E cambian entre fases: los módulos se emparejaron por **identidad celular mapeada dentro de ese día**. En SAMPLE–TEST comparten 24/33 y 20/40 células core (`J=0,727` y `0,500`; `q≈0,013` ambos). Veinte y diecinueve células core, respectivamente, aparecen en las cuatro fases.

   - Módulo A: **OF1 E1 → SAMPLE E2 → TEST E1 → OF2 E1**.
   - Módulo B: **OF1 E2 → SAMPLE E1 → TEST E2 → OF2 E3**.

Si hay espacio, usá estas figuras como **exploración metodológica o apoyo**. Una matriz binaria de coseno sirve para explicar el método, pero no es la conclusión biológica; el timeline del ajuste global es descriptivo porque los factores se definieron usando las cuatro fases a la vez.

## Si te preguntan «¿qué es un ensemble y qué es el core?»

- **Ensemble:** un patrón recurrente de vectores de actividad simultánea; Stoixeion lo obtiene de una matriz de similitud coseno y SVD. Tiene momentos de activación y células miembro.
- **Core neurons:** el subconjunto de células más representativo de ese patrón. Una célula puede participar en más de un ensemble. E1 de SAMPLE y E1 de TEST son etiquetas locales de dos ajustes independientes; se comparan **sus células**, no sus números.
- **Qué prueba la estadística:** que el solapamiento de **membresía core** entre SAMPLE y TEST supera el de identidades celulares permutadas. No es una prueba independiente de que cada episodio temporal se haya reproducido, de memoria del objeto S ni de tratamiento.
- **Por qué el `q` viene del ZIP anterior:** la corrida Fast reciente usa 99 permutaciones auxiliares (`p` mínimo 0,01) y sus pares no pasan FDR. La corrida phasewise anterior usó 999; para las nueve jornadas focales, sus 2.589 membresías core y 42.774 asignaciones de ensemble a frame coinciden exactamente con Fast. Usamos **Fast para la visualización nueva** y **999 para el valor `q`**, sin mezclarlos.

## Si preguntan por CNO o por el Factor 5 de R005

«Vemos pares core relacionados en jornadas VEH y CNO, pero **no medimos un efecto farmacológico uniforme**. R004 CNO es la excepción del resumen S–D. El Factor global 5 de R005 VEH es llamativo porque casi todas sus asignaciones caen en TEST; lo mostraría solo como una pista individual. El día CNO también tiene actividad tardía de TEST. El reporte focal del cluster respalda por huella `A` **2/18** enlaces core entre fechas; las otras 16 son indeterminadas. Su mapping global recuperó solo 83/136 células 4/4 VEH y 97/165 CNO, así que no sirve para comparar ensembles poblacionales».

## Frase y leyenda en inglés para el póster

**Exploratory method result:** “Core cell sets of CA1 coactivity patterns identified independently in SAMPLE and TEST overlapped above a shuffled-identity null in five of six S–D sessions after multiple-comparison correction. This supports within-day recurrence of ensemble membership; it does not establish object-specific coding, temporal replay, or a treatment effect.”

**Caption:** “Each point shows the strongest SAMPLE–TEST core pair per session. Jaccard overlap was assessed against 999 cell-identity permutations; Benjamini–Hochberg correction covered 176 candidate SAMPLE–TEST pairs across ten object sessions. Filled symbols pass `q≤0.05`; the open symbol does not. Six sessions are nested in three mice.”

## Separación con el otro análisis

El resultado de **coactividad por pares de neuronas** (6/6 S–D frente a dos nulos) es un análisis distinto. No digas que ese `p` valida los ensembles de Stoixeion. Si te preguntan específicamente por Carrillo-Reid, mostrá el solapamiento de cores y el ejemplo temporal de R006. Las listas completas de factores y episodios permanecen en el paquete local curado.
