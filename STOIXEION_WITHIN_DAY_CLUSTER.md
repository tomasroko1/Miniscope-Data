# Stoixeion dentro de cada fecha

Este flujo corre en MATLAB sobre los `.mat` diarios de `DataBase/<ANIMAL>`. Usa el `act.mapping` que ya está dentro de cada archivo para alinear células entre sus subsesiones. Corre Stoixeion por separado en cada fase y omite el ajuste concatenado de cuatro fases y cualquier registro entre fechas.

En MATLAB, con el repo ya actualizado, pegá esta línea en la ventana de comandos:

```matlab
addpath('/mnt/NAS/Tomas/Miniscope-Data'); run_stoixeion_miniscope('phasewise','/mnt/NAS/Miniscopes/Reg_CA1/DataBase','/mnt/NAS/Tomas/Stoixeion','/mnt/NAS/Tomas/results/stoixeion/phasewise')
```

El selector `phasewise` incluye HabL y las jornadas SD/XsS de cuatro fases. Cada fase genera sus figuras y resultados; `core_overlap.csv` resume el solapamiento de los cores entre fases. `E1` de una fase no se considera automáticamente el mismo ensemble que `E1` de otra: el análisis compara los cores mediante las identidades celulares que aporta `act.mapping`.

Cada carpeta de fase incluye `singular_values.png`: espectro de la SVD de la matriz de similitud, con los rangos seleccionados marcados en rojo. Los valores num?ricos tambi?n quedan en `singular_values.csv` en la carpeta ra?z de resultados.

Los `.mat` simples de HabC1/HabC2 no contienen un mapping de cuatro fases y no entran en esta corrida.
