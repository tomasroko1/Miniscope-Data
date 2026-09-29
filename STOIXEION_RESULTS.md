# Resultados de Stoixeion

Ubicacion central: `/mnt/NAS/Tomas/results/stoixeion/`. El enlace `current`
apunta a la corrida mas reciente completada. Las carpetas `runs/` son corridas
independientes: no se deben sumar sus CSV porque tienen alcances y parametros
distintos.

| Carpeta | Fecha | Alcance | Estado |
| --- | --- | --- | --- |
| `current` -> `runs/2026-09-28_fast_9_focal_sessions` | 28/9/2026 | R004, R005 y R006; HabL, SD CNO y SD VEH; por fase y concatenado dentro de cada fecha; 99 shuffles | 41 resultados `ok` de 44 registros. Tres fases R004 HabL no superaron el criterio de umbral. |
| `phasewise` -> `runs/2026-09-27_phasewise_14_sessions` | 27/9/2026 | 14 sesiones, incluye XsS; solo por fase | 49 fases `ok` de 56. Tres fases R004 HabL sin umbral significativo y cuatro fases R004 XsS VEH con error de `RasNum`. |
| `archive/` | Anterior | Pilotos, intento fallido y analisis global con mappings dudosos | Historico; no usar como resultado vigente. |

Ambas corridas usan el `act.mapping` interno de cada `.mat` para las cuatro
fases de una fecha. `GLOBAL_CONCAT` en la corrida nueva significa concatenar
esas cuatro fases dentro de la misma fecha; no usa los mappings globales
externos entre fechas.

## Donde encontrar cada cosa

En `current/workers/R004`, `R005` y `R006` estan los CSV originales de cada
worker y sus logs. `current/tables/` contiene los CSV consolidados. Los CSV en
la raiz de `current/` son enlaces a `tables/` para mantener las rutas previas.
`current/figures/` contiene un panel resumen, ocho mapas de recurrencia y
ocho graficos de actividad global, regenerados de los CSV sin repetir el
analisis. Esta corrida desactivo los PNG por fase durante MATLAB.

La corrida `phasewise` anterior conserva sus figuras por fase, incluidos los
paneles de actividad poblacional y los espectros singulares.

Las rutas viejas `/mnt/NAS/Tomas/results/stoixeion_fast_clean` y
`Miniscope-Data/results/stoixeion/phasewise` siguen funcionando como enlaces.
En el repositorio estan los runners `run_fast_cluster.sh`,
`run_stoixeion_fast.m` y `run_stoixeion_miniscope.m`; la dependencia Stoixeion
se prepara con `setup_stoixeion.sh` y `stoixeion_exports.patch`.
