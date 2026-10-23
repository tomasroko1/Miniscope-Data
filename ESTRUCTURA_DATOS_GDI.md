# Documentación y Diccionario de Datos: Pipeline `GDi/` y Miniscope CA1

Este documento detalla la jerarquía completa de archivos, scripts y estructuras de datos en el directorio de trabajo del laboratorio:
`/mnt/NAS/Miniscopes/Reg_CA1/GDi`

Abarca desde el directorio raíz (`GDi/`), pasando por las carpetas de los animales (`R004_M4_26`, `R005_M5_26`, `R006_M6_26`), las sesiones experimentales (`Hab C1`, `Hab C2`, `Hab L`, `T1` a `T4`), los bloques de tarea (`A_AV1`, `B_SAMPLE`, `C_TEST`, `D_AV2`), hasta las carpetas de resultados detallados como `CaliAli_results/`, `My_V4_Miniscope/` y `My_WebCam/`.

---

## 1. Nivel Raíz: `GDi/`

Ruta: `/mnt/NAS/Miniscopes/Reg_CA1/GDi`

Contiene la colección de experimentos con miniscopios inalámbricos/alámbricos (UCLA V4 Miniscope) y cámaras de comportamiento en la región CA1 del hipocampo, junto con las herramientas centrales de procesamiento cinemático, control de calidad morfométrico y tracking por DeepLabCut (DLC).

### Carpetas de Sujetos / Animales
- **`R001`, `R002`, `R003`, `R007`, `R008`**: Directorios de cohortes de animales adicionales.
- **`R004_M4_26`**: Cohorte 2026, Ratón 4, Miniscope 4.
- **`R005_M5_26`**: Cohorte 2026, Ratón 5, Miniscope 5.
- **`R006_M6_26`**: Cohorte 2026, Ratón 6, Miniscope 6.

### Scripts de Procesamiento y Control de Calidad en `GDi/`
- **`analizar_animales.py`**:
  Orquestador principal de tracking de comportamiento. Rastrea recursivamente las carpetas de los animales en busca de videos en `My_WebCam/`, comprueba si existen inferencias previas de DeepLabCut (`.h5`, `.csv`) y lanza la ejecución sobre los videos pendientes.
- **`analizar_animales_shuffle2.py`** / **`analizar_animales_cambiar suffle 1.py`**:
  Variantes de `analizar_animales.py` configuradas para trabajar con diferentes conjuntos de entrenamiento (*shuffles* 1 o 2 de DeepLabCut) y modelos de redes neuronales (ej. HRNet-W32 vs ResNet).
- **`analisis_geometria.py`** y **`score_geometria.py`**:
  Módulos de validación morfofisiológica del tracking. Calculan las distancias euclidianas entre pares anatómicos adyacentes (`Nose-Neck`, `Neck-Spine`, `Spine-TailB`, `TailB-TailE`). Generan métricas de distorsión corporal (*score de geometría*) para detectar y marcar automáticamente saltos espaciales espurios o fallos de inferencia de DLC.
- **`filter_tracking_v1.py`** / **`filter_tracking_v2_diagnostic.py`** / **`filter_tracking_mark.py`** / **`mark_v2_diagnostic.py`**:
  Pipelines de filtrado y limpieza de trayectorias. Aplican filtros de mediana temporal, interpolación de fotogramas con baja confianza de probabilidad ($p < \text{umbral}$) y restricciones biomecánicas sobre las coordenadas corporales.
- **`prueba_filtro_temporal.py`**:
  Script de prueba unitaria y ajuste de hiperparámetros de suavizado temporal sobre series de tiempo de posición $(x, y)$.
- **`clasificar_eventos_v2.py`**:
  Clasificador de conducta animal basado en cinemática (cálculo de velocidad instantánea, aceleración, vector de orientación cefálica e interacción/aproximación con objetos).
- **`inspeccionar_eventos_v2.py`**:
  Herramienta diagnóstica para auditar e inspeccionar gráficamente los eventos conductuales clasificados.
- **`visualizar_marcas.py`**:
  Script de renderizado interactivo que superpone en tiempo real sobre el video de la cámara cenital el esqueleto predicho por DeepLabCut, las marcas filtradas y los límites de la arena y objetos.
- **`Using`**:
  Archivo marcador temporal de estado de ejecución (0 bytes).
- **`__pycache__/`**:
  Archivos compilados de Python (`.pyc`).

---

## 2. Nivel Animal: `R004_M4_26`, `R005_M5_26`, `R006_M6_26`

Representan los datos longitudinales completos de cada animal individual.

### Convención de Nombres
- **`R00X`**: Identificador único del animal / registro.
- **`MX`**: Dispositivo Miniscope utilizado (ej. M4, M5, M6).
- **`26`**: Año de adquisición experimental (2026).

### Scripts y Logs del Pipeline de Imágenes de Calcio
- **`run_all_days_minian.py`** / **`run_all_days_minian_v2.py`**:
  Script de automatización por animal que itera sobre todas las carpetas de sesiones del sujeto y lanza de manera ordenada el procesamiento por lotes de MiniAn.
- **`run_all_days_minian.log`**:
  Archivo de registro completo con los tiempos de cómputo, uso de memoria, advertencias y convergencia de las corridas de MiniAn en todos los días.
- **`run_minian_batch.py`**:
  Pipeline principal de MiniAn: realiza la corrección de movimiento del sensor (*rigid motion correction*), remoción de fondo neuropilar y extracción de fuentes neuronales (CNMF) para aislar las trazas de fluorescencia $\Delta F/F$ y eventos discretos.
- **`run_cross_registration.py`** / **`run_cross_registration_v2.py`**:
  Algoritmo de registro cruzado inter-bloques o inter-días (*cross-registration*). Realiza la alineación espacial rígida/afín de los centroides neuronales y calcula las matrices de coincidencia celular para seguir a la misma neurona a lo largo del tiempo.
- **`06_export_mappings_and_C_for_matlab.py`** *(en R004)*:
  Script para exportar las matrices de actividad $C$ de MiniAn junto con los índices de correspondencia celular a estructuras compatibles con MATLAB (`.mat`), utilizadas en el paquete `Stoixeion`.
- **`cross_registration_log.log`** / **`minian_batch.log`**:
  Logs específicos de las etapas de registro celular y extracción por lotes.

---

## 3. Nivel Sesión: Carpetas Experimentales

Cada carpeta dentro del animal representa un día o hito experimental:

### Tipos de Sesión
1. **Habituación a Campo Abierto**:
   - `Hab C1_2026_06_16` / `Hab C2_2026_06_17`:
     Fases de familiarización al entorno de campo abierto circular/cuadrado sin objetos estimulares.
2. **Habituación a Pista Lineal**:
   - `Hab L_2026_06_18` *(o `Hab L_2026_06_18_stripes_in_11_and_12`)*:
     Familiarización a desplazamiento unidimensional / pista lineal con marcas visuales de referencia.
3. **Fases de Test Conductual / Farmacogenética**:
   - `T1_SD_...` a `T4_XsS_...`:
     - **`SD`** (*Same Day*): Paradigma donde la muestra (Sample) y la prueba (Test) se realizan en el mismo día con intervalo corto.
     - **`XsS`** (*Across Sessions / Cross-Session*): Paradigma de retención a largo plazo con separación temporal entre fases.
     - **`CNO`**: Administración intraperitoneal previa de Clozapina-N-óxido (activador de receptores sintéticos DREADDs para modulación/silenciamiento neuronal específico de CA1).
     - **`VEH`**: Administración de vehículo / solución salina estéril como condición control intra-sujeto.

### Archivos de la Sesión
- **`cents.pkl`**: Coordenadas espaciales de los centroides somáticos de las neuronas identificadas en esa jornada.
- **`mappings.csv`** / **`mappings.pkl`**: Matriz de mapeo que asigna qué neurona de cada bloque (`A`, `B`, `C`, `D`) corresponde a la misma identidad celular global.
- **`mappings_matlab.csv`** / **`mappings_for_matlab.csv`**: Versión con indexación 1-based para MATLAB.
- **`shiftds.nc`**: Campo de desplazamiento bidimensional (formato NetCDF / xarray) derivado de la alineación para corregir derivas del miniscope.
- **`cross_registration_alignment.png`**: Imagen de diagnóstico que superpone con código de colores los mapas celulares de cada bloque.
- **`cell_contours.csv.gz`**: Coordenadas comprimidas de los polígonos que delimitan los contornos de cada ROI neuronal.
- **`cross_registration_contours.pdf`**: Reporte gráfico vectorial de los contornos celulares sobre el campo de visión.
- **`unit_id_to_C_column.csv`**: Tabla de traducción que vincula el identificador global de la célula (`unit_id`) con su columna correspondiente en la matriz de fluorescencia $C$.
- **`open_field_alignment.mat`**: Calibración geométrica y límites de la arena espacial exportados para análisis cinemático en MATLAB.
- **`run_*_USED_*.py`**: Snapshots congelados del código exacto ejecutado en esa sesión histórica para trazabilidad y reproducibilidad.
- **`05_temporal_sparse_registration_sweep.py`** / **`temporal_sparse_registration_sweep.log`**: Barridos de evaluación de la densidad del registro neuronal a lo largo del tiempo.

---

## 4. Nivel Bloque: Subsesiones (`A_AV1`, `B_SAMPLE`, `C_TEST`, `D_AV2`)

En los días de test, cada sesión se subdivide en 4 fases temporales:
- **`A_AV1`**: Ambiente Vacío 1 (registro basal previo al test).
- **`B_SAMPLE`**: Fase de Muestra (exploración de objetos novedosos).
- **`C_TEST`**: Fase de Test (evaluación de discriminación/memoria).
- **`D_AV2`**: Ambiente Vacío 2 (registro posterior para control de deriva).

### Estructura Interna de Cada Bloque:
1. **`My_V4_Miniscope/`**:
   - `0.avi`, `1.avi`, ...: Videos secuenciales crudos en escala de grises de la señal de fluorescencia en CA1.
   - `timeStamps.csv`: Estampas de tiempo (en ms) de cada fotograma capturado por el sensor CMOS del miniscopio.
   - `headOrientation.csv`: Trazas temporales de orientación cefálica (roll, pitch, yaw, cuaterniones) del sensor IMU integrado.
   - `minian/`: Carpeta Zarr con matrices intermedias de preprocesamiento, corrección de movimiento y CNMF.
   - `minian_C_for_matlab.mat`: Matriz de fluorescencia de calcio denoised ($C$) exportada a MATLAB.
2. **`My_WebCam/`**:
   - `0.avi`: Grabación de video de alta definición de la cámara cenital del comportamiento del animal.
   - `*DLC*.h5`, `*DLC*.csv`, `*DLC*.pickle`: Predicciones de coordenadas $(x, y, p)$ para cada punto anatómico y objeto generadas por DeepLabCut.
3. **Metadatos y Hardware**:
   - `metaData.json`: Configuración del miniscopio (ganancia, exposición, potencia del LED de excitación).
   - `notes.csv`: Anotaciones del operador durante la captura.
   - `behaviorTracker/` y `experiment/`: Configuraciones del software de adquisición de UCLA.

---

## 5. Nivel Especializado: `CaliAli_results/`

Ubicado dentro de cada sesión experimental en `R005_M5_26` y `R006_M6_26`.

**CaliAli** (*Calcium Alignment*) es una suite desarrollada en MATLAB diseñada para realizar el registro inter-sesión elástico de campos neuronales, la extracción no negativa de señales (CNMF) y la cuantificación morfofisiológica de eventos de calcio.

### Contenido de `CaliAli_results/`

#### 1. `CaliAli_traces_<ANIMAL>_<SESION>.mat`
Formato MATLAB v7.3 (HDF5) que contiene las series temporales neuronales:
- **`traces.C`**: Matriz de trazas de fluorescencia desconvueltas y filtradas ($\Delta F/F$ limpias).
- **`traces.C_raw`**: Trazas de fluorescencia crudas sin corrección de deconvolución.
- **`traces.S`**: Raster de picos o eventos discretos de calcio inferidos (*deconvolved spike raster*).
- **`traces.session`**: Vector que clasifica a qué bloque (`A_AV1`, `B_SAMPLE`, `C_TEST`, `D_AV2`) pertenece cada muestra temporal.
- **`traces.raw_frame`**: Índices originales de los fotogramas en los videos crudos.
- **`traces.stripe_frames`**: Vector de fotogramas donde se detectaron artefactos de franjas de ruido (*stripes*).
- **`fps`**: Frecuencia de muestreo del registro (Hz).
- **`traces_note`**: Metadatos explicativos sobre el algoritmo de cálculo de trazas.

#### 2. `CaliAli_summary_<ANIMAL>_<SESION>.csv` y `.mat`
Tabla resumen a nivel de neurona individual. Cada fila corresponde a una célula unificada:
- **`neuron`**: Identificador de la neurona.
- **`centroid_x`, `centroid_y`**: Coordenadas espaciales $(X, Y)$ del centroide del soma.
- **`area_px`**: Área del soma en píxeles.
- **`n_sessions_active`**: Cantidad de sub-bloques en los que la neurona estuvo activa.
- **Variables por fase** (`A_AV1`, `B_SAMPLE`, `C_TEST`, `D_AV2`):
  - `*_events`: Conteo total de eventos de calcio registrados en esa fase.
  - `*_events_min`: Tasa de disparo / frecuencia de eventos por minuto.
  - `*_snr`: Relación Señal-Ruido (*Signal-to-Noise Ratio*) calculada para esa célula.
  - `*_active`: Flag binario (1 o 0) que define si la célula estuvo activa según el umbral mínimo de eventos.
- **Variables adicionales en `.mat`**:
  - `contours_xy`: Polígonos vectoriales con las coordenadas del perímetro de cada soma.
  - `contour_level`: Nivel umbral utilizado para la segmentación del contorno.
  - `correlation_image`: Matriz bidimensional de la correlación espacial local de píxeles.
  - `dims_hw`: Alto y ancho del sensor en píxeles.

#### 3. `CaliAli_options_<ANIMAL>_<SESION>.mat`
Contiene la estructura `CaliAli_options` con los parámetros aplicados:
- **`downsampling`**: Factores de submuestreo espacial y temporal.
- **`preprocessing`**: Filtros pasa-banda espaciales y sustracción morfológica de fondo.
- **`motion_correction`**: Parámetros del algoritmo de corrección de movimiento intra-bloque (NoRMCorre).
- **`inter_session_alignment`**: Configuración de los algoritmos de optimización de traslación y registro afín entre bloques.
- **`cnmf`**: Hiperparámetros de factorización CNMF (radio esperado de neuronas $gSig$, umbral de aceptación, regularización).

#### 4. `CaliAli_result_<ANIMAL>_<SESION>.mat`
Archivo binario maestro (entre 200 MB y 850 MB) con la estructura completa `neuron`:
- Matrices de huella espacial $A$ (distribución de intensidad en el sensor para cada célula).
- Matrices de fondo espacial $b$ y fondo temporal $f$ (neuropilo difuso).
- Residuos de la deconvolución y convergencia matemática del modelo.

#### 5. `CaliAli_alignment_<ANIMAL>_<SESION>.txt`
Reporte cuantitativo en texto plano con las métricas de registro inter-bloque:
- **`BV score`** (*Background Variation score*): Métrica de estabilidad del fondo entre sesiones.
- **Tabla comparativa de alineación**: Compara tres estados:
  1. *Non-aligned*: Sin corrección inter-bloque.
  2. *Translations*: Corrección de desplazamientos rígidos en $X$ e $Y$.
  3. *CaliAli*: Alineamiento elástico no rígido completo.
- **Métricas de salida**:
  - `Correlation scores`: Correlación espacial entre las proyecciones alineadas de cada bloque.
  - `Mean Corr. Score`: Correlación promedio global del campo celular.
  - `Crispness`: Nitidez/contraste espacial de la imagen promedio alineada (a mayor valor, menor desenfoque y mejor superposición).
  *(En sesiones de un solo bloque como `Hab C1` y `Hab C2`, reporta que no aplica alineamiento inter-sesión).*

#### 6. Diagnósticos Visuales (`*.png`)
- **`CaliAli_alignment_*_projections.png`**: Imagen diagnóstica de las proyecciones compuestas del campo de visión antes y después de la alineación para inspeccionar visualmente distorsiones o rotaciones.
- **`CaliAli_alignment_*_scores.png`**: Curvas de evolución de las métricas de correlación y optimización espacial.
- **`CaliAli_summary_*.png`**: Panel principal de control de calidad con la imagen de correlación del tejido, contornos celulares superpuestos e histogramas de SNR, actividad y eventos.
