# Stoixeion en el cluster

Desde la raíz de `Miniscope-Data`, actualizá el repo y prepará Stoixeion una sola vez:

```bash
git pull origin main
if [ ! -d ../Stoixeion ]; then
  git clone https://github.com/luiscareid/Stoixeion.git ../Stoixeion
fi
if ! grep -q "Diagnostics = struct();" ../Stoixeion/Stoixeion.m; then
  git -C ../Stoixeion apply --unidiff-zero --ignore-space-change "$PWD/stoixeion_exports.patch"
fi
export MINISCOPE_DATA_DIR="$PWD/data"
export STOIXEION_DIR="$PWD/../Stoixeion"
```

Probá primero una jornada. Para sacar rápido el diccionario global de todas las jornadas, corré `global`. Si queda tiempo, corré `phasewise` como control sin repetir el análisis global; `primary` hace ambos en una sola corrida:

```bash
matlab -batch "addpath(pwd); run_stoixeion_miniscope('pilot')"
matlab -batch "addpath(pwd); run_stoixeion_miniscope('global')"
matlab -batch "addpath(pwd); run_stoixeion_miniscope('phasewise')"
```

`MINISCOPE_DATA_DIR` debe ser la carpeta que contiene `R004`, `R005` y `R006`; cambiá esa ruta si tus datos están en otro lugar.

`global` procesa HabL, SD y XsS solo con el análisis concatenado. `phasewise` hace solo los cuatro análisis separados. `primary` procesa esas mismas jornadas con ambos análisis. Para limitarlo: `habl`, `sd` o `xss`. `all4` intenta todas las sesiones de cuatro fases, con ambos análisis.

En cada jornada el runner hace dos análisis: uno independiente por fase y otro con las cuatro fases concatenadas. Para el global usa las neuronas mapeadas en las cuatro fases y un único umbral por neurona, `S > 3 × SD(S)` calculado sobre las cuatro fases; las corridas fasewise calculan su umbral en cada fase. Stoixeion detecta los factores globales una vez y luego el CSV informa cuántos vectores significativos de cada factor aparecen por fase. Los números de factor solo identifican factores dentro de esa jornada; no enlazan días.

Resultados en `results/stoixeion/<seleccion>/`: CSV agregados de fases, cores, miembros, valores singulares, umbrales, activaciones y actividad de factores globales por fase. Las figuras quedan por animal/jornada; dentro de `GLOBAL_CONCATENATED` está `global_factor_phase_activity.png`. No se generan `.npz`.

La permutación nativa de Stoixeion en la corrida concatenada mezcla tiempos de las cuatro fases y no conserva las tasas específicas de cada fase. Además, fases con más vectores pueden pesar más al construir el diccionario global; el CSV informa tasas por minuto para compararlas. Por eso el global es exploratorio: contrastalo con los análisis independientes por fase y no interpretes una diferencia como efecto CNO/VEH sin resumir por animal.
