# Seguimiento para el agente del cluster

Ya hiciste el panorama inicial y creaste `results/CA1_poster_curado_2026_09_29/`. **Continuá desde ese trabajo; no repitas el inventario.** Este handoff aporta el runner de registro y el export trazable de coactividad S–D. Respondé en español. Conservá las salidas anteriores.

## Prioridad hoy: cerrar las afirmaciones del póster

1. Verificá `network_6day/PRIMARY_6DAY_VERIFICATION.csv` contra los CSV por jornada, protocolos y `audit_sd_pairwise_network.py` incluidos en este paquete. Son seis jornadas S–D; con onsets `S>3 SD`, bins de 1 s y modo global retirado, las seis superan nulos de identidad y corrimiento (`p=0,005` por jornada; Holm `p=0,03`). `python audit_sd_pairwise_network.py --self-test` pasó en la PC local. Copiá al índice curado la tabla, figura y procedencia exacta. **Esta comprobación no incluye control de posición/velocidad.**
2. El resultado **10/10 con control de posición/velocidad** está aquí solo como tabla agregada y figura (`network_10day_reported/`). El código `make_ca1_context_network_figure.py` cita cuatro CSV detallados bajo `All_object_networks/` que **faltan en el repositorio local**. Buscalos en el cluster. Si aparecen, verificá sus parámetros, 10 filas, nulos, p ajustadas y hashes, y adjuntalos al índice. Si no aparecen, etiquetá el 10/10 como *reportado, procedencia detallada pendiente* y usá el 6/6 verificable para el póster.
3. Presentá R005 Factor 5 como **pista intradía**: 959 frames en TEST bajo VEH, 2 en SAMPLE y 2 en OF2; en CNO hay actividad tardía en TEST incorporada al Factor 2 global. El QC de huella `A` reportado hasta ahora es **0/18 células core concluyentes**. Dos enlaces tabulares y seis candidatos por centroide no prueban conservación ni desaparición del ensemble. No uses la comparación entre días como resultado farmacológico del póster.

## Reanudar el registro focal de R005

El intento anterior usó `/mnt/NAS/Tomas/run_cross_registration_by_animal.py`; esa ruta es incorrecta para el checkout documentado. Comprobá primero `/mnt/NAS/Tomas/Miniscope-Data/run_cross_registration_by_animal.py`. Si no existe o es otra versión, este handoff incluye el runner actualizado y sus dos módulos auxiliares en `registration_runner/`; usalos juntos desde una carpeta de trabajo nueva. El manifiesto de las ocho adquisiciones también está ahí.

El intento anterior también falló por ausencia de `xarray`. Antes de lanzar MiniAn, identificá un entorno que importe **minian, xarray, dask, skimage, scipy, pandas y matplotlib**. Revisá los entornos conda existentes, incluido `/mnt/NAS/Tomas/miniconda311`, y usá uno que pase una prueba de importación. Si ninguno sirve, creá un entorno nuevo siguiendo `R004_REGISTRATION_CLUSTER.md` del repo y comprobá las importaciones. Evitá instalar paquetes sueltos en un entorno de otro análisis.

Con runner y entorno listos, corré desde una carpeta nueva y única, con las rutas reales del cluster:

```bash
python /mnt/NAS/Tomas/Miniscope-Data/run_cross_registration_by_animal.py --animal-root /mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26 --session-list /mnt/NAS/Tomas/Miniscope-Data/registration_manifests/R005_SD_VEH_CNO.txt --reference-mat-dir /mnt/NAS/Miniscopes/Reg_CA1/DataBase/R005 --output-dir /mnt/NAS/Tomas/results/R005_SD_VEH_CNO_registration_targeted_20260929_v2 --require-qc-figures
```

Si el script o manifiesto no están en el checkout, apuntá el comando a las **copias del handoff**; no improvises otra selección de sesiones. Revisá que las 136 células VEH y 165 CNO 4/4 mantengan sus pares intradía. Inspeccioná contornos `A` de todos los enlaces posibles de las 18 células core F5 y exigí una lista con `matched / rechazado / indeterminado`, distancia de centroide, solapamiento de huella y motivo. Indicá el número de células realmente evaluables; no trates las no enlazadas como desaparecidas. Si el runner no produce el QC de contornos, documentá esa limitación y no asciendas candidatos por centroide a identidades confirmadas.

## Actualización que necesito

Actualizá `INDEX.md`, `R005_CROSSDAY.md` y `PANTALLAZO_PARA_REVISAR.md` en la carpeta curada. Tu respuesta final debe decir, en este orden: **(a)** qué afirmación de coactividad puede ir hoy al póster con fuente auditable; **(b)** si el registro focal terminó y cuántas células core F5 tienen identidad entre días confirmada por huella; **(c)** qué queda fuera del póster; **(d)** rutas de los nuevos CSV, figuras y logs. Si el registro sigue bloqueado, dame el error exacto y la siguiente acción concreta; terminá de todos modos la curaduría del 6/6.
