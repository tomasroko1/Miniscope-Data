# Flujo CaliAli v5

El usuario identificó `/mnt/NAS/Miniscopes/CaliAli/v5/` como el flujo que usa actualmente. `v5/` contiene una copia de sus archivos `caliali_*.m` y de su README original, obtenida el 13 de noviembre de 2026. El archivo `CaliAli_my_parameters.m` se copió de la carpeta padre `/mnt/NAS/Miniscopes/CaliAli/`: es una referencia; comprobar cuál está en el MATLAB path de cada corrida. El ejemplo original es `v5/caliali_workflow_example.m`.

La copia **no incluye el paquete CaliAli-main** ni sus dependencias. El servidor usa una instalación externa con parches locales; el README original enumera los parches. No está probado que esta copia sea ejecutable fuera del entorno original. No asumir que v5 generó todos los archivos históricos.

## Etapas según el ejemplo original

1. `caliali_status(animal)`: consultar el estado de cada día.
2. `caliali_run(animal, 'mc')`: downsampling, tratamiento de franjas y corrección de movimiento.
3. `caliali_crop_animal(animal)`: **interactivo**, un recorte compartido para el animal.
4. `caliali_run(animal, 'align')`: alineación de las fases; revisar sus imágenes de control.
5. `caliali_set_seeds(day)`: **interactivo**, selección de umbrales y máscara por día.
6. `caliali_run(animal, 'extract')`: extracción CNMF-E y resumen para los días preparados.

Los resultados exportados incluyen `CaliAli_result`, `CaliAli_summary` y `CaliAli_traces`. En las trazas, las filas son neuronas y las columnas cuadros; las fases se guardan separadas. Se incluyen `C`, `C_raw`, `S`, `raw_frame` y las marcas de cuadros con franjas. `raw_frame` es una referencia importante para sincronizar con comportamiento, pero esa sincronización aún necesita validación.

## Efectos sobre archivos

Los comandos originales trabajan sobre rutas configuradas del NAS y del directorio de trabajo. No son operaciones de solo lectura: según la etapa, crean, reemplazan o eliminan archivos intermedios. En particular, `caliali_clean_day` elimina archivos y `caliali_redeconv` con un único valor puede reescribir resúmenes y trazas compartidas. **No ejecutar estas operaciones sobre los originales para una prueba aislada.**

Para pruebas de redeconvolución, preparar una salida aislada y verificar las trazas antes de sustituir cualquier resultado compartido.
