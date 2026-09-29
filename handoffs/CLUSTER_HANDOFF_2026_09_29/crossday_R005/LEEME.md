# R005: registro entre VEH y CNO, 19 y 21 de junio

## Resultado corto

Encontré en Descargas el ZIP `R005-20260927T013505Z-1-001.zip`, exportado de un registro MiniAn de **22 adquisiciones** de R005 que terminó en el cluster el 26/9. Audité sus dos mappings contra `act.mapping` de los `.mat` diarios y ejecuté un emparejamiento **focal de los dos días S–D** usando los centroides alineados exportados. Copias del [mapping corregido](source_exports/mappings_preserving_existing.csv), el [mapping original](source_exports/mappings_22sessions_original.csv), los [centroides](source_exports/cents_aligned.csv) y el [log](source_exports/cross_registration_log.log) están separadas de las tablas derivadas.

| Registro | Consistencia con las células presentes en las 4 fases VEH/CNO | Células enlazadas entre ambos días 4/4 | Factor 5 VEH con enlace a CNO |
| --- | ---: | ---: | ---: |
| `mappings.csv` de las 22 adquisiciones | **18/136 y 18/165**; muchos conflictos dentro de cada día | 9 | No utilizable para esta pregunta |
| `mappings_preserving_existing.csv` del ZIP | **136/136 y 164/165** | **9** | **2/18** células core; ambas forman parte del core del Factor 2 CNO |
| Nuevos candidatos geométricos focales, solo centroides | Parte de los 136/165 IDs diarios 4/4 | **46 candidatos**, incluidos los 9 enlaces anteriores | **6/18** candidatas; 4 están en el core Factor 2 CNO |

**Interpretación:** los dos enlaces ya presentes y consistentes por tablas son compatibles con que una parte del core de Factor 5 VEH participe en Factor 2 CNO. Son **candidatos de identidad**, todavía sin huella `A` confirmada. Las otras 46 parejas son **candidatas de posición**; no permiten afirmar identidad del ensemble ni cuantificar su conservación. Por tanto, el registro disponible no demuestra que Factor 5 desaparezca bajo CNO. La actividad tardía de TEST observada en el ajuste CNO sigue siendo una explicación alternativa importante.

**Actualización del cluster, 29/9 (reporte focal más reciente):** terminó un registro de las **ocho adquisiciones** VEH/CNO. Su mapping global **no es utilizable para conclusiones poblacionales**: conserva 83/136 células VEH y 97/165 CNO presentes en las cuatro fases dentro de cada día. El QC local de huellas `A` respalda **2/18** enlaces de células core del Factor 5 entre fechas; ambas aparecen en el core del Factor 2 CNO y su ensemble de TEST. Las otras 16 siguen **indeterminadas**, no ausentes. El intento anterior que falló por ruta/`xarray` queda superado. La auditoría focal y sus tablas están únicamente en el cluster, bajo `/mnt/NAS/Tomas/results/R005_SD_VEH_CNO_registration_targeted_20260929_v2/`; aquí solo dispongo del resumen mostrado por el agente, no de sus CSV para reauditarlo.

## Cómo se produjo la lista de 46 candidatas

Cada una de las 136 células VEH y 165 CNO de `act.mapping` presentes en las cuatro fases recibió el **centroide promedio** de sus cuatro adquisiciones, usando los centroides ya alineados del ZIP. Una pareja entra solo si los centroides de las dos fechas están a **≤3 píxeles**, son vecinos más próximos **mutuos** y cada segundo vecino está al menos **2 píxeles** más lejos en ambos días. El criterio es uno-a-uno. Los nueve enlaces 4/4 del mapping corregido están dentro de esas 46 parejas. [Imagen de QC](R005_crossday_centroid_qc.png), [parejas candidatas](centroid_candidates_3px_margin2.csv), [solo las seis candidatas del Factor 5](F5_core_centroid_candidates.csv), [sensibilidad 2–5 px](centroid_candidate_sensitivity.csv).

La mediana de la desviación máxima de los cuatro centroides respecto del promedio dentro de un día es 2,58 px en VEH y 2,41 px en CNO. Los nueve enlaces existentes tienen distancias entre centros de 0,33–1,99 px. El alineamiento general se ve razonable en la proyección del ZIP, pero **un centroide cercano no prueba que dos huellas sean la misma neurona**, sobre todo en CA1 denso. El ZIP no contiene los contornos ni matrices `A` para hacer esa validación local.

## Las dos células core enlazadas en el registro existente

| ID de célula dentro del día VEH | ID dentro del día CNO | ID del mapping entre días | Distancia entre centroides |
| ---: | ---: | ---: | ---: |
| 14 | 362 | 3 | 1,72 px |
| 300 | 36 | 828 | 1,24 px |

Ambas están en el core del Factor 5 VEH y del Factor 2 CNO según las tablas de mapping. El denominador importante es **2 de 18** células core con enlace 4/4 **tabular**; el reporte focal posterior también respalda **2 de 18** por huella `A`. Sin los CSV del cluster aquí no puedo afirmar independientemente que sean exactamente las mismas dos filas de esta tabla. Las otras 16 no se pueden tratar como ausentes del día CNO. [Los nueve enlaces](existing_consistent_all4_links.csv) y [el subconjunto de dos](F5_core_existing_links.csv) están en CSV. Los `mapped_cell_id` diarios empiezan en **0**; el `global_cell_id` del nuevo mapping empieza en **1**.

## Registro focal: comando archivado y resultado posterior

El ajuste de 22 adquisiciones no es la prueba óptima para estos dos días. El cluster tiene los `A` originales y **ya ejecutó** el registro focal de ocho adquisiciones VEH/CNO, con contornos y QC. El manifiesto exacto derivado de `sess.sess_paths` está en [R005_SD_VEH_CNO.txt](R005_SD_VEH_CNO.txt). La ruta de trabajo que muestran los scripts es:

```text
Repositorio: /mnt/NAS/Tomas/Miniscope-Data
MiniAn:      /mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26
MAT:         /mnt/NAS/Miniscopes/Reg_CA1/DataBase/R005
```

El siguiente comando queda archivado como instrucción de reproducción; **no hace falta repetirlo para el póster**:

```bash
python /mnt/NAS/Tomas/Miniscope-Data/run_cross_registration_by_animal.py \
  --animal-root /mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26 \
  --session-list /mnt/NAS/Tomas/Miniscope-Data/registration_manifests/R005_SD_VEH_CNO.txt \
  --reference-mat-dir /mnt/NAS/Miniscopes/Reg_CA1/DataBase/R005 \
  --output-dir /mnt/NAS/Tomas/results/R005_SD_VEH_CNO_registration_targeted_20260929 \
  --require-qc-figures
```

El control posterior encontró que solo 83/136 y 97/165 IDs diarios 4/4 se conservaron en el mapping global. **No usar** ese `mappings.csv` para comparar poblaciones o factores entre días. El QC focal de las 18 células core es una pregunta separada y solo respalda dos enlaces individuales.

La fuente exacta y hashes SHA-256 están en [source_sha256.json](source_sha256.json). Esta PC dispone de los `.mat` y del ZIP exportado, pero no tiene `/mnt/NAS` montado ni SSH configurado; por eso aquí no es posible repetir la etapa de registro de huellas `A`.
