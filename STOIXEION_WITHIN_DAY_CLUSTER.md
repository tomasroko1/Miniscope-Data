# Stoixeion por fase dentro de cada fecha

Este es el flujo usado para la corrida phasewise de R004, R005 y R006. Lee los
`.mat` de `DataBase/<ANIMAL>`, usa `act.mapping` dentro de cada archivo para
identificar las celulas comunes entre las cuatro fases de esa fecha, y ejecuta
Stoixeion por separado en cada fase. No alinea celulas entre fechas ni usa los
mappings globales externos.

## Codigo y dependencias

- `run_stoixeion_miniscope.m` prepara los datos, llama a Stoixeion y exporta
  tablas y figuras. La seleccion `phasewise` incluye HabL y las jornadas SD/XsS.
- `stoixeion_exports.patch` contiene todos los cambios hechos al codigo
  original de Stoixeion para esta corrida. El original es
  `luiscareid/Stoixeion` en el commit
  `be6f547cfb6d2b49c1c0ec3b6331d9ea98f08491`.
- `setup_stoixeion.sh` descarga ese commit y aplica el parche. Por defecto
  prepara `external/Stoixeion` dentro de este repo; tambien acepta otra ruta
  como primer argumento. No sobrescribe cambios locales en un clon existente.
- Se uso MATLAB R2017a con Statistics and Machine Learning Toolbox (`pdist2`
  y `perfcurve`).

Preparar la dependencia y ejecutar el lote completo:

```bash
bash setup_stoixeion.sh
/usr/local/MATLAB/R2017a/bin/matlab -nodisplay -nodesktop -nosplash -r "try, addpath(pwd); run_stoixeion_miniscope('phasewise','/mnt/NAS/Miniscopes/Reg_CA1/DataBase',fullfile(pwd,'external','Stoixeion'),'/mnt/NAS/Tomas/results/stoixeion/phasewise'); catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"
```

Ejecutar el comando desde la raiz de `Miniscope-Data`. Elegir un directorio de
salida nuevo para una nueva corrida; el runner escribe CSV incrementales y
archivos de figura dentro de ese directorio.

La corrida de septiembre de 2026 se dividio por animal con
`run_stoixeion_animal_worker.sh` y se reunio en
`/mnt/NAS/Tomas/results/stoixeion/phasewise` con
`consolidate_stoixeion_phasewise.py`. Ese script de consolidacion corresponde
a la estructura de carpetas especifica de aquella corrida; no es necesario
para un lote nuevo ejecutado directamente con el comando anterior.

## Resultados y figuras

`phase_summary.csv` registra el estado de cada fase. `event_thresholds.csv`
guarda los umbrales de binarizacion por celula; `singular_values.csv` y
`ensemble_activity.csv` guardan el espectro y las apariciones de ensambles.
Los identificadores de ensamble son locales a cada fase: E1 en dos fases
distintas no implica el mismo ensamble. `core_overlap.csv` compara los cores
entre fases usando `act.mapping`.

Cada fase completada tiene `stoixeion_04.png` con dos paneles: actividad
poblacional y raster de apariciones por ensamble. Los numeros de ensamble ya
no se dibujan sobre la actividad. `singular_values.png` muestra solo los
rangos que Stoixeion considera para detectar ensambles, con eje Y logaritmico.

Si los CSV y MAT ya existen, se pueden regenerar estas figuras sin repetir el
analisis, el shuffling ni la SVD:

```bash
python3 regenerate_stoixeion_population_activity.py
/usr/local/MATLAB/R2017a/bin/matlab -nodisplay -nodesktop -nosplash -r "try, addpath(pwd); plot_stoixeion_singular_values('/mnt/NAS/Tomas/results/stoixeion/phasewise'); catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"
```

La regeneracion de actividad usa NumPy, SciPy y Matplotlib. Ambos scripts de
regeneracion tienen rutas por defecto de la corrida en este NAS y se deben
adaptar si los datos se mueven.
