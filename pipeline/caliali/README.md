# Procesamiento con CaliAli v5

El punto de entrada es [`v5/caliali_workflow_example.m`](v5/caliali_workflow_example.m), una copia del ejemplo que se usa en el clúster. Se ejecuta **por secciones** desde MATLAB; el recorte y la selección de semillas son interactivos.

## Dependencias

El archivo de ejecución usa las funciones y la instalación de CaliAli que ya están en el clúster. No las incluye dentro del repositorio:

```text
/mnt/NAS/Miniscopes/CaliAli/v5/       Funciones del flujo
/mnt/data1/CaliAli/CaliAli-main/      Instalación CaliAli del servidor
```

La configuración `CaliAli_my_parameters.m` debe estar disponible en el entorno de MATLAB. El ejemplo del clúster usa un directorio de trabajo local y guarda los resultados finales en la carpeta de cada sesión en el NAS. Comprobar rutas y parámetros en el servidor donde se ejecute. Esta copia del repositorio no fue validada como instalación independiente.

## Orden de trabajo

1. **Estado:** elegir el animal y consultar sus días con `caliali_status`.
2. **Corrección:** ejecutar `caliali_run(animal, 'mc')` para preparar los videos y corregir movimiento.
3. **Recorte:** ejecutar `caliali_crop_animal(animal)` y seleccionar el área interactiva.
4. **Alineación:** ejecutar `caliali_run(animal, 'align')` y revisar las imágenes de control.
5. **Semillas:** ejecutar `caliali_set_seeds` para cada día y ajustar sus parámetros en la interfaz.
6. **Extracción:** ejecutar `caliali_run(animal, 'extract')` para los días preparados.

Los resultados incluyen `CaliAli_result`, `CaliAli_summary` y `CaliAli_traces`. Las trazas contienen `C`, `C_raw`, `S`, los nombres de fase y referencias a cuadros originales y artefactos.

## Antes de ejecutar

El script de ejemplo apunta a rutas compartidas y sus etapas pueden crear, reemplazar o borrar archivos intermedios. La limpieza y la aplicación de un valor individual de redeconvolución también pueden modificar resultados. Para pruebas, usar una copia aislada y una salida nueva. Conservar el mismo conjunto de parámetros entre días del animal cuando se planee registrarlos entre sí.
