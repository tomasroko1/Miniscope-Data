# Seis días S–D: resultado del atlas local

Se reconstruyeron **25 factores globales** de **seis jornadas** (tres animales × VEH/CNO). Los seis `.mat` pasaron el QC de mapping, umbral, número de vectores significativos y frames asignados. El ajuste SVD de los factores proviene del export MATLAB existente; las curvas, matrices muestreadas y tablas se calcularon aquí desde los `.mat`.

| Animal | Tratamiento | Células 4/4 | Factores | Core <3 células | Frames asignados |
| --- | --- | ---: | ---: | ---: | ---: |
| R004 | VEH | 99 | 6 | 3 | 3297 |
| R004 | CNO | 169 | 6 | 6 | 2113 |
| R005 | VEH | 136 | 5 | 0 | 4497 |
| R005 | CNO | 165 | 4 | 1 | 4967 |
| R006 | VEH | 65 | 2 | 0 | 5292 |
| R006 | CNO | 83 | 2 | 0 | 8037 |

**Control de interpretabilidad:** 24/25 factores incluyen como «miembros» a más del 80 % de las células mapeadas; 10/25 tienen un core de una o dos células. Esas banderas son descriptivas, no un filtro estadístico. La lista *core* resulta más informativa que la membresía amplia para comparar grupos celulares. En R004 CNO, los seis factores tienen core de una o dos células, por lo que no los mostraría como ensambles de muchas neuronas sincronizadas.

**Ejemplo destacado, R005 VEH Factor 5:** 0/2/959/2 frames asignados en OF1/SAMPLE/TEST/OF2. Es un patrón de ese día concentrado en TEST. El día CNO también presenta actividad tardía en TEST, pero sus factores se ajustaron de nuevo y no comparten identidad garantizada con VEH.

**Conclusión:** el atlas resuelve qué factores se detectaron y cuándo se expresan dentro de cada jornada. No identifica qué información espacial representan ni demuestra un efecto CNO. Para la metodología, mirar el panel de coseno, espectro y asignaciones; para la dinámica, las figuras altas de cada día.
