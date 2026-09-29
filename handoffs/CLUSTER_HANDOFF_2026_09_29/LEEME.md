# Handoff al agente del cluster, 29/9/2026

Pasale al agente el archivo `PROMPT_REANUDAR_CLUSTER.md` de esta carpeta y el ZIP que la contiene. El agente ya devolvió un panorama inicial; este paquete sirve para **resolver sus dos faltantes** sin repetir ese inventario.

- `registration_runner/`: versión local de `run_cross_registration_by_animal.py`, sus dos módulos auxiliares y el manifiesto exacto de ocho adquisiciones R005 VEH/CNO. En el intento anterior se llamó una ruta equivocada del script y faltó `xarray` en el entorno.
- `network_6day/`: código, protocolos, CSV por jornada, figura y [resumen primario verificable](network_6day/PRIMARY_6DAY_VERIFICATION.csv) de las seis jornadas S–D. Las seis superan los dos nulos tras retirar el modo global; son tres animales. La prueba sintética del script pasó localmente. Esta versión **no controla posición/velocidad**.
- `network_10day_reported/`: tabla agregada y figura de diez jornadas, con el script que las ensambló. En esta PC faltan los cuatro CSV detallados `All_object_networks/` que el script necesita; hasta recuperarlos, el 10/10 con posición/velocidad queda **pendiente de auditoría de procedencia**.
- `crossday_R005/`: informe y tablas locales de enlaces tabulares y candidatos por centroides. El panorama del cluster informó **0/18 células core F5 con QC concluyente de huella A**. Ninguna tabla por centroides sustituye esa comprobación.
- `carrillo_poster/`: dos figuras Stoixeion listas para explicar el solapamiento de cores de las seis jornadas S–D y la dinámica de dos módulos de R006 CNO. Su `q` proviene de la corrida phasewise de **999 permutaciones**; es un análisis distinto de la huella de red por pares.
- `SHA256_MANIFEST.csv`: ruta de origen, tamaño y hash de cada archivo copiado. No se modificaron los originales.

El póster puede apoyarse ya en el **6/6 S–D trazable** y mostrar R005 Factor 5 como ejemplo intradía. La validación interdiaria del Factor 5 y el resultado controlado 10/10 quedan fuera de las afirmaciones principales hasta completar sus controles.
